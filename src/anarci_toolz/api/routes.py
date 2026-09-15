from fastapi import APIRouter, Depends, HTTPException
from starlette.concurrency import run_in_threadpool

from .dependencies import get_semaphore
from .models import NumberRequest, NumberResponse
from .settings import Settings, get_settings

router = APIRouter()


@router.post("/number", response_model=NumberResponse)
async def number(
    request: NumberRequest,
    settings: Settings = Depends(get_settings),
) -> NumberResponse:
    if len(request.sequences) > settings.MAX_SEQUENCES:
        raise HTTPException(
            status_code=400,
            detail=(
                f"Received {len(request.sequences)} sequences, which exceeds "
                f"the maximum of {settings.MAX_SEQUENCES} per request."
            ),
        )

    semaphore = get_semaphore()

    try:
        async with semaphore:
            # TODO(jon): wire up real call to anarci_toolz.pipeline.run_anarci_toolz
            # once the seq_dna -> DataFrame translation is figured out.
            # run_anarci_toolz is CPU-bound/multiprocessing, so it must run via
            # run_in_threadpool to avoid blocking the event loop, e.g.:
            #   df_result = await run_in_threadpool(run_anarci_toolz, df, ...)
            results: list[dict] = await run_in_threadpool(lambda: [])
            return NumberResponse(results=results)
    except HTTPException:
        raise
    except Exception:
        raise HTTPException(
            status_code=500,
            detail="An unexpected error occurred while processing the request.",
        )
