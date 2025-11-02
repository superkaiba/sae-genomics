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

from sae_lens.training.sae import SAE
from sae_lens.training.train_sae_on_activations import train_sae_on_activations

from sae_genomics.models.tahoe_adapter import TahoeModelAdapter


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

    def initialize_sae(self):
        """Initialize the SAE model."""
        self.sae = SAE(
            d_in=self.d_in,
            d_sae=self.d_sae,
            activation_fn=self.activation,
            k=self.k if self.activation == "topk" else None,
        )
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
        metrics = {"losses": [], "l1_losses": [], "recon_losses": []}

        for epoch in range(n_epochs):
            pbar = tqdm(dataloader, desc=f"Epoch {epoch+1}/{n_epochs}")

            for batch_idx, (batch_acts,) in enumerate(pbar):
                batch_acts = batch_acts.to(self.device)

                # Forward pass
                sae_out, feature_acts, loss_data = self.sae(
                    batch_acts,
                    return_loss=True,
                    l1_coefficient=self.l1_coefficient,
                )

                loss = loss_data["loss"]
                l1_loss = loss_data.get("l1_loss", 0)
                recon_loss = loss_data.get("reconstruction_loss", 0)

                # Backward pass
                optimizer.zero_grad()
                loss.backward()
                optimizer.step()

                # Logging
                if total_steps % log_every == 0:
                    metrics["losses"].append(loss.item())
                    metrics["l1_losses"].append(l1_loss if isinstance(l1_loss, float) else l1_loss.item())
                    metrics["recon_losses"].append(recon_loss if isinstance(recon_loss, float) else recon_loss.item())

                    pbar.set_postfix({
                        "loss": f"{loss.item():.4f}",
                        "l1": f"{metrics['l1_losses'][-1]:.4f}",
                        "recon": f"{metrics['recon_losses'][-1]:.4f}",
                    })

                # Checkpointing
                if output_dir is not None and total_steps % save_every == 0 and total_steps > 0:
                    self.save_checkpoint(output_dir / f"checkpoint_step_{total_steps}.pt")

                total_steps += 1

        # Save final checkpoint
        if output_dir is not None:
            self.save_checkpoint(output_dir / "final_checkpoint.pt")

        return metrics

    def train_on_dataloader(
        self,
        tahoe_adapter: TahoeModelAdapter,
        dataloader: DataLoader,
        n_steps: int = 10000,
        batch_size: int = 256,
        extract_every: int = 1000,
        log_every: int = 100,
        output_dir: Optional[Path] = None,
    ) -> Dict:
        """Train SAE by extracting activations on-the-fly.

        Args:
            tahoe_adapter: Tahoe model adapter
            dataloader: DataLoader with single-cell data
            n_steps: Total training steps
            batch_size: SAE training batch size
            extract_every: Extract new activations every N steps
            log_every: Log frequency
            output_dir: Directory to save checkpoints

        Returns:
            Training metrics dictionary
        """
        if self.sae is None:
            self.initialize_sae()

        # Extract initial activations
        print("Extracting initial activations...")
        activations, _ = tahoe_adapter.extract_activations(
            dataloader,
            max_batches=extract_every // dataloader.batch_size,
        )

        # Flatten
        activations_flat = activations.reshape(-1, activations.shape[-1])

        # Create dataset
        dataset = TensorDataset(activations_flat)
        sae_dataloader = DataLoader(dataset, batch_size=batch_size, shuffle=True)

        # Setup optimizer
        optimizer = torch.optim.Adam(self.sae.parameters(), lr=self.lr)

        # Training loop
        self.sae.train()
        total_steps = 0
        metrics = {"losses": [], "l1_losses": [], "recon_losses": []}

        pbar = tqdm(total=n_steps, desc="Training SAE")

        while total_steps < n_steps:
            # Re-extract activations periodically
            if total_steps % extract_every == 0 and total_steps > 0:
                print(f"\nRe-extracting activations at step {total_steps}...")
                activations, _ = tahoe_adapter.extract_activations(
                    dataloader,
                    max_batches=extract_every // dataloader.batch_size,
                )
                activations_flat = activations.reshape(-1, activations.shape[-1])
                dataset = TensorDataset(activations_flat)
                sae_dataloader = DataLoader(dataset, batch_size=batch_size, shuffle=True)

            # Train on current activation buffer
            for batch_acts, in sae_dataloader:
                if total_steps >= n_steps:
                    break

                batch_acts = batch_acts.to(self.device)

                # Forward pass
                sae_out, feature_acts, loss_data = self.sae(
                    batch_acts,
                    return_loss=True,
                    l1_coefficient=self.l1_coefficient,
                )

                loss = loss_data["loss"]

                # Backward pass
                optimizer.zero_grad()
                loss.backward()
                optimizer.step()

                # Logging
                if total_steps % log_every == 0:
                    metrics["losses"].append(loss.item())
                    pbar.set_postfix({"loss": f"{loss.item():.4f}"})

                # Checkpointing
                if output_dir is not None and total_steps % 1000 == 0 and total_steps > 0:
                    self.save_checkpoint(output_dir / f"checkpoint_step_{total_steps}.pt")

                total_steps += 1
                pbar.update(1)

        pbar.close()

        # Save final checkpoint
        if output_dir is not None:
            self.save_checkpoint(output_dir / "final_checkpoint.pt")

        return metrics

    def save_checkpoint(self, path: Path):
        """Save SAE checkpoint."""
        path.parent.mkdir(parents=True, exist_ok=True)
        torch.save({
            "sae_state_dict": self.sae.state_dict(),
            "d_in": self.d_in,
            "d_sae": self.d_sae,
            "activation": self.activation,
            "k": self.k,
        }, path)
        print(f"Saved checkpoint to {path}")

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

        trainer.initialize_sae()
        trainer.sae.load_state_dict(checkpoint["sae_state_dict"])

        print(f"Loaded checkpoint from {path}")
        return trainer
