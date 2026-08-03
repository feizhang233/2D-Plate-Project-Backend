"""FastAPI wrapper for the Kirchhoff/Mindlin plate mathematics."""

from __future__ import annotations

import json
import os
from importlib.resources import files
from pathlib import Path
from uuid import uuid4

import uvicorn
from fastapi import FastAPI, HTTPException, Request
from fastapi.responses import JSONResponse
from fastapi.staticfiles import StaticFiles

from .api_models import PlateAnalysisRequest, PlateAnalysisResponse
from .plotting import render_postprocess_images
from .service import PlateInputError, serialize_analysis, solve_plate_case


def _template() -> dict:
    resource = files("mindlin_plate").joinpath("templates/rectangular_plate.json")
    return json.loads(resource.read_text(encoding="utf-8"))


def create_app(output_root: str | Path | None = None) -> FastAPI:
    """Application factory; ``output_root`` is injectable for tests/deployment."""

    root = Path(
        output_root
        or os.environ.get("MINDLIN_PLATE_OUTPUT_DIR", "outputs/plate-analyses")
    ).resolve()
    root.mkdir(parents=True, exist_ok=True)
    application = FastAPI(
        title="2D Plate Analysis API",
        version="1.0.0",
        description=(
            "FastAPI service for automatic Kirchhoff DKQ / Reissner-Mindlin "
            "MITC4 plate analysis and PNG post-processing."
        ),
    )
    application.state.output_root = root
    application.state.results = {}
    application.mount("/artifacts", StaticFiles(directory=str(root)), name="artifacts")

    @application.exception_handler(PlateInputError)
    async def plate_input_error(_: Request, exc: PlateInputError) -> JSONResponse:
        return JSONResponse(status_code=422, content={"detail": str(exc)})

    @application.exception_handler(RuntimeError)
    async def postprocess_error(_: Request, exc: RuntimeError) -> JSONResponse:
        return JSONResponse(status_code=500, content={"detail": str(exc)})

    @application.get("/health", tags=["system"])
    async def health() -> dict[str, str]:
        return {"status": "ok"}

    @application.get("/api/v1/templates/rectangular-plate", tags=["preprocess"])
    async def rectangular_plate_template() -> dict:
        """Return a ready-to-edit preprocessing JSON document."""

        return _template()

    @application.post(
        "/api/v1/analyses",
        response_model=PlateAnalysisResponse,
        status_code=201,
        tags=["analysis"],
    )
    def create_analysis(payload: PlateAnalysisRequest) -> dict:
        case = payload.model_dump(mode="json")
        analysis = solve_plate_case(case)
        analysis_id = uuid4().hex
        analysis_directory = root / analysis_id
        paths = render_postprocess_images(
            analysis,
            analysis_directory,
            plots=payload.postprocess.plots,
            dpi=payload.postprocess.dpi,
        )
        response = serialize_analysis(analysis)
        response.update(
            {
                "analysis_id": analysis_id,
                "status": "completed",
                "images": [
                    {
                        "field": path.stem,
                        "media_type": "image/png",
                        "url": f"/artifacts/{analysis_id}/{path.name}",
                    }
                    for path in paths
                ],
            }
        )
        application.state.results[analysis_id] = response
        return response

    @application.post(
        "/api/v1/solve",
        response_model=PlateAnalysisResponse,
        status_code=201,
        include_in_schema=False,
    )
    def solve_alias(payload: PlateAnalysisRequest) -> dict:
        """Short alias retained for command-line and integration clients."""

        return create_analysis(payload)

    @application.get(
        "/api/v1/analyses/{analysis_id}",
        response_model=PlateAnalysisResponse,
        tags=["analysis"],
    )
    async def get_analysis(analysis_id: str) -> dict:
        response = application.state.results.get(analysis_id)
        if response is None:
            raise HTTPException(status_code=404, detail="analysis not found")
        return response

    return application


app = create_app()


def main() -> None:
    host = os.environ.get("MINDLIN_PLATE_HOST", "127.0.0.1")
    port = int(os.environ.get("MINDLIN_PLATE_PORT", "8000"))
    uvicorn.run("mindlin_plate.api:app", host=host, port=port)


if __name__ == "__main__":
    main()
