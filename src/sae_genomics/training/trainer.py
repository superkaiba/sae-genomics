"""SAE training on Tahoe X1 activations."""

import sys
from pathlib import Path
from typing import Dict, Optional

import torch
from torch.utils.data import DataLoader, TensorDataset
from tqdm import tqdm

# Add SAELens to path
SAELENS_PATH = Path(__file__).parent.parent.parent.parent / "external" / "SAELens"
if str(SAELENS_PATH) not in sys.path:
    sys.path.insert(0, str(SAELENS_PATH))

from sae_lens import TopKTrainingSAE, TopKTrainingSAEConfig, TopKSAE, TopKSAEConfig

from sae_genomics.models.tahoe_adapter import TahoeModelAdapter

# Training constants
GRADIENT_CLIP_VALUE = 1.0  # Max gradient norm for clipping


class SAETrainer:
    """Trainer for sparse autoencoders on Tahoe X1 activations."""

    def __init__(
        self,
        d_in: int,
        d_sae: int,
        activation: str = "topk",
        k: int = 64,
        lr: float = 3e-4,
        l1_coefficient: float = 1e-3,
        device: str = "auto",
    ):
        """Initialize SAE trainer.

        Args:
            d_in: Input dimension (model d_model)
            d_sae: SAE dictionary size
            activation: Activation function ('topk' or 'relu')
            k: Number of active features for topk
            lr: Learning rate
            l1_coefficient: L1 sparsity penalty
            device: Device to train on
        """
        self.d_in = d_in
        self.d_sae = d_sae
        self.activation = activation
        self.k = k
        self.lr = lr
        self.l1_coefficient = l1_coefficient

        if device == "auto":
            self.device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
        else:
            self.device = torch.device(device)

        # Initialize SAE
        self.sae = None
        self.sae_inference = None  # For inference after training

    def initialize_sae(self):
        """Initialize the SAE model."""
        config = TopKTrainingSAEConfig(
            d_in=self.d_in,
            d_sae=self.d_sae,
            k=self.k,
            device=str(self.device),
        )
        self.sae = TopKTrainingSAE(config)
        self.sae = self.sae.to(self.device)

    def train_on_activations(
        self,
        activations: torch.Tensor,
        n_epochs: int = 10,
        batch_size: int = 256,
        log_every: int = 100,
        save_every: int = 1000,
        output_dir: Optional[Path] = None,
    ) -> Dict:
        """Train SAE on pre-extracted activations.

        Args:
            activations: Tensor of shape (n_samples, seq_len, d_model)
            n_epochs: Number of training epochs
            batch_size: Batch size
            log_every: Log frequency
            save_every: Checkpoint save frequency
            output_dir: Directory to save checkpoints

        Returns:
            Training metrics dictionary
        """
        if self.sae is None:
            self.initialize_sae()

        # Flatten activations: (n_samples * seq_len, d_model)
        n_samples, seq_len, d_model = activations.shape
        activations_flat = activations.reshape(-1, d_model)

        # Create dataset and dataloader
        dataset = TensorDataset(activations_flat)
        dataloader = DataLoader(
            dataset,
            batch_size=batch_size,
            shuffle=True,
            num_workers=0,  # Use 0 for tensor dataset
        )

        # Setup optimizer
        optimizer = torch.optim.Adam(self.sae.parameters(), lr=self.lr)

        # Training loop
        self.sae.train()
        total_steps = 0
        metrics = {"losses": [], "mse_losses": [], "aux_losses": []}

        # Calculate total training steps for SAELens
        total_training_steps = len(dataloader) * n_epochs

        # Coefficients for training (L1 coefficient, etc.)
        coefficients = {
            "l1": self.l1_coefficient,
        }

        for epoch in range(n_epochs):
            pbar = tqdm(dataloader, desc=f"Epoch {epoch+1}/{n_epochs}")

            for batch_idx, (batch_acts,) in enumerate(pbar):
                batch_acts = batch_acts.to(self.device)

                # Create TrainStepInput with required parameters
                from sae_lens.saes.sae import TrainStepInput

                step_input = TrainStepInput(
                    sae_in=batch_acts,
                    coefficients=coefficients,
                    dead_neuron_mask=None,  # Simplified for now
                    n_training_steps=total_training_steps,
                )

                # Forward pass
                output = self.sae.training_forward_pass(step_input)

                # Total loss is in output.losses dictionary
                total_loss = output.losses.get("mse_loss", torch.tensor(0.0))
                for loss_name, loss_value in output.losses.items():
                    if loss_name != "mse_loss" and isinstance(loss_value, torch.Tensor):
                        total_loss = total_loss + loss_value

                # Backward pass
                optimizer.zero_grad()
                total_loss.backward()

                # Gradient clipping for training stability
                torch.nn.utils.clip_grad_norm_(self.sae.parameters(), GRADIENT_CLIP_VALUE)

                optimizer.step()

                # Logging
                if total_steps % log_every == 0:
                    mse = output.losses.get("mse_loss", torch.tensor(0.0)).item()
                    aux = sum(
                        v.item()
                        for k, v in output.losses.items()
                        if k != "mse_loss" and isinstance(v, torch.Tensor)
                    )

                    metrics["losses"].append(total_loss.item())
                    metrics["mse_losses"].append(mse)
                    metrics["aux_losses"].append(aux)

                    pbar.set_postfix({
                        "loss": f"{total_loss.item():.4f}",
                        "mse": f"{mse:.4f}",
                    })

                # Checkpointing
                if output_dir is not None and total_steps % save_every == 0 and total_steps > 0:
                    self.save_checkpoint(output_dir / f"checkpoint_step_{total_steps}.pt")

                total_steps += 1

        # Save final checkpoint
        if output_dir is not None:
            self.save_checkpoint(output_dir / "final_checkpoint.pt")

        return metrics

    def save_checkpoint(self, path: Path):
        """Save SAE checkpoint."""
        path.parent.mkdir(parents=True, exist_ok=True)

        # Convert to inference SAE
        config = TopKSAEConfig(
            d_in=self.d_in,
            d_sae=self.d_sae,
            k=self.k,
            device=str(self.device),
        )
        sae_inference = TopKSAE(config)

        # Copy weights - filter to only include keys that exist in inference SAE
        # (Training SAE has auxiliary loss parameters that inference SAE doesn't need)
        state_dict = self.sae.state_dict()
        inference_keys = set(sae_inference.state_dict().keys())
        filtered_state_dict = {k: v for k, v in state_dict.items() if k in inference_keys}

        # Load with strict=False and check for issues
        missing_keys, unexpected_keys = sae_inference.load_state_dict(filtered_state_dict, strict=False)

        # Missing keys in inference SAE are a critical error
        if missing_keys:
            raise RuntimeError(
                f"Failed to convert training SAE to inference SAE. "
                f"Missing required keys: {missing_keys}"
            )

        # Unexpected keys are fine (they're training-only parameters)
        if unexpected_keys:
            print(f"Note: Skipped training-only parameters: {unexpected_keys}", file=sys.stderr)

        torch.save({
            "sae_state_dict": sae_inference.state_dict(),
            "d_in": self.d_in,
            "d_sae": self.d_sae,
            "activation": self.activation,
            "k": self.k,
        }, path)
        print(f"Saved checkpoint to {path}", file=sys.stderr)

    @classmethod
    def load_checkpoint(cls, path: Path, device: str = "auto") -> "SAETrainer":
        """Load SAE from checkpoint."""
        checkpoint = torch.load(path, map_location="cpu")

        trainer = cls(
            d_in=checkpoint["d_in"],
            d_sae=checkpoint["d_sae"],
            activation=checkpoint["activation"],
            k=checkpoint.get("k", 64),
            device=device,
        )

        # Load as inference SAE
        config = TopKSAEConfig(
            d_in=checkpoint["d_in"],
            d_sae=checkpoint["d_sae"],
            k=checkpoint.get("k", 64),
            device=str(trainer.device),
        )
        trainer.sae_inference = TopKSAE(config)
        trainer.sae_inference.load_state_dict(checkpoint["sae_state_dict"])
        trainer.sae_inference = trainer.sae_inference.to(trainer.device)

        # Use inference SAE as main SAE
        trainer.sae = trainer.sae_inference

        import sys
        print(f"Loaded checkpoint from {path}", file=sys.stderr)
        return trainer
