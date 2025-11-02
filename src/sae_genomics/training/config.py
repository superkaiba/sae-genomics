"""Training configuration utilities."""

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Optional

import yaml


@dataclass
class TrainingConfig:
    """Configuration for SAE training."""

    # Experiment
    experiment_name: str = "sae_training"
    run_name: Optional[str] = None

    # Model
    model_config: str = "configs/models/tx1_70m.yaml"

    # SAE architecture
    d_in: int = 512
    d_sae: int = 16384
    expansion_factor: int = 32
    activation: str = "topk"
    k: int = 64

    # Training
    batch_size: int = 32
    total_training_steps: int = 10000
    lr: float = 3e-4
    l1_coefficient: float = 1e-3

    # Logging
    log_every: int = 100
    save_every: int = 1000
    eval_every: int = 500

    # Data
    data_path: Optional[str] = None
    max_cells: Optional[int] = None

    # Device
    device: str = "auto"
    seed: int = 42

    @classmethod
    def from_yaml(cls, path: Path) -> "TrainingConfig":
        """Load configuration from YAML file."""
        with open(path) as f:
            config = yaml.safe_load(f)

        training_cfg = config.get("training", {})
        sae_cfg = training_cfg.get("sae", {})
        optimizer_cfg = training_cfg.get("optimizer", {})
        data_cfg = training_cfg.get("data", {})

        return cls(
            experiment_name=training_cfg.get("experiment_name", "sae_training"),
            run_name=training_cfg.get("run_name"),
            model_config=training_cfg.get("model_config", "configs/models/tx1_70m.yaml"),
            d_in=sae_cfg.get("d_in", 512),
            d_sae=sae_cfg.get("d_sae", 16384),
            expansion_factor=sae_cfg.get("expansion_factor", 32),
            activation=sae_cfg.get("activation", "topk"),
            k=sae_cfg.get("k", 64),
            batch_size=training_cfg.get("batch_size", 32),
            total_training_steps=training_cfg.get("total_training_steps", 10000),
            lr=optimizer_cfg.get("lr", 3e-4),
            l1_coefficient=training_cfg.get("loss", {}).get("l1_coefficient", 1e-3),
            log_every=training_cfg.get("log_every", 100),
            save_every=training_cfg.get("save_every", 1000),
            eval_every=training_cfg.get("eval_every", 500),
            data_path=data_cfg.get("dataset_path"),
            max_cells=data_cfg.get("max_cells"),
            device=training_cfg.get("device", "auto"),
            seed=training_cfg.get("seed", 42),
        )

    def to_dict(self) -> Dict:
        """Convert to dictionary."""
        return {
            "experiment_name": self.experiment_name,
            "run_name": self.run_name,
            "model_config": self.model_config,
            "d_in": self.d_in,
            "d_sae": self.d_sae,
            "expansion_factor": self.expansion_factor,
            "activation": self.activation,
            "k": self.k,
            "batch_size": self.batch_size,
            "total_training_steps": self.total_training_steps,
            "lr": self.lr,
            "l1_coefficient": self.l1_coefficient,
            "log_every": self.log_every,
            "save_every": self.save_every,
            "eval_every": self.eval_every,
            "data_path": self.data_path,
            "max_cells": self.max_cells,
            "device": self.device,
            "seed": self.seed,
        }
