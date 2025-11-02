"""Model configuration utilities."""

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional

import yaml


@dataclass
class ModelConfig:
    """Configuration for Tahoe X1 models."""

    name: str
    size: str
    n_params: int
    context_length: int
    d_model: int
    n_layers: int
    n_heads: int
    checkpoint_path: Optional[str] = None
    from_pretrained: bool = True
    cache_dir: str = "./data/models"
    device: str = "auto"
    dtype: str = "float32"
    hook_points: Optional[List[str]] = None

    @classmethod
    def from_yaml(cls, path: Path) -> "ModelConfig":
        """Load configuration from YAML file."""
        with open(path) as f:
            config = yaml.safe_load(f)

        model_cfg = config.get("model", {})
        return cls(
            name=model_cfg.get("name", "tahoe-x1"),
            size=model_cfg.get("size", "70M"),
            n_params=model_cfg.get("n_params", 70_000_000),
            context_length=model_cfg.get("context_length", 1024),
            d_model=model_cfg.get("d_model", 512),
            n_layers=model_cfg.get("n_layers", 12),
            n_heads=model_cfg.get("n_heads", 8),
            checkpoint_path=model_cfg.get("checkpoint_path"),
            from_pretrained=model_cfg.get("from_pretrained", True),
            cache_dir=model_cfg.get("cache_dir", "./data/models"),
            device=model_cfg.get("device", "auto"),
            dtype=model_cfg.get("dtype", "float32"),
            hook_points=model_cfg.get("hook_points"),
        )

    def to_dict(self) -> Dict:
        """Convert to dictionary."""
        return {
            "name": self.name,
            "size": self.size,
            "n_params": self.n_params,
            "context_length": self.context_length,
            "d_model": self.d_model,
            "n_layers": self.n_layers,
            "n_heads": self.n_heads,
            "checkpoint_path": self.checkpoint_path,
            "from_pretrained": self.from_pretrained,
            "cache_dir": self.cache_dir,
            "device": self.device,
            "dtype": self.dtype,
            "hook_points": self.hook_points,
        }

    def get_hf_model_size(self) -> str:
        """Get HuggingFace model size string."""
        if "70m" in self.size.lower():
            return "70m"
        elif "1.3b" in self.size.lower() or "1b" in self.size.lower():
            return "1.3b"
        elif "3b" in self.size.lower():
            return "3b"
        else:
            raise ValueError(f"Unknown model size: {self.size}")
