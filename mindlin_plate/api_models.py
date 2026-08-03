"""Pydantic transport models for the FastAPI application."""

from __future__ import annotations

from typing import Literal

from pydantic import BaseModel, ConfigDict, Field


class StrictModel(BaseModel):
    model_config = ConfigDict(extra="forbid")


class UnitsInput(StrictModel):
    length: str = "m"
    force: str = "N"
    stress: str = "Pa"


class MeshInput(StrictModel):
    type: Literal["rectangular", "custom"] = "rectangular"
    length_x: float | None = Field(default=None, gt=0.0)
    length_y: float | None = Field(default=None, gt=0.0)
    elements_x: int | None = Field(default=None, ge=1, le=200)
    elements_y: int | None = Field(default=None, ge=1, le=200)
    nodes: list[list[float]] | None = None
    elements: list[list[int]] | None = None


class MaterialInput(StrictModel):
    young: float = Field(gt=0.0)
    poisson: float = Field(gt=-1.0, lt=0.5)
    thickness: float = Field(gt=0.0)
    shear_correction: float = Field(default=5.0 / 6.0, gt=0.0)


class LoadInput(StrictModel):
    type: Literal["uniform", "sinusoidal"] = "uniform"
    magnitude: float


class BoundaryInput(StrictModel):
    type: Literal[
        "clamped",
        "soft_simply_supported",
        "hard_simply_supported",
        "custom",
    ] = "hard_simply_supported"
    prescribed_dofs: dict[int, float] | None = None


class AnalysisOptions(StrictModel):
    plate_method: Literal["auto", "K", "M"] = "auto"
    shear_scheme: Literal["full", "reduced", "mitc4"] = "mitc4"
    thinness_threshold: float = Field(default=1.0 / 20.0, gt=0.0, lt=1.0)


class PostprocessOptions(StrictModel):
    plots: list[Literal["deflection", "rotation", "moment", "stress_top"]] = Field(
        default_factory=lambda: ["deflection", "rotation", "moment", "stress_top"]
    )
    dpi: int = Field(default=140, ge=72, le=300)


class PlateAnalysisRequest(StrictModel):
    name: str = Field(default="plate-analysis", min_length=1, max_length=120)
    units: UnitsInput = Field(default_factory=UnitsInput)
    mesh: MeshInput
    material: MaterialInput
    load: LoadInput
    boundary: BoundaryInput = Field(default_factory=BoundaryInput)
    analysis: AnalysisOptions = Field(default_factory=AnalysisOptions)
    postprocess: PostprocessOptions = Field(default_factory=PostprocessOptions)


class ImageArtifact(BaseModel):
    field: str
    media_type: Literal["image/png"] = "image/png"
    url: str


class PlateAnalysisResponse(BaseModel):
    analysis_id: str
    status: Literal["completed"] = "completed"
    name: str
    summary: dict
    mesh: dict
    nodal_results: list[dict]
    element_results: list[dict]
    images: list[ImageArtifact]
