import asyncio
from typing import Optional

from .settings import get_settings

# run_parallel_anarci's Pool(processes=cpu_count()) assumes a single job owns
# the whole node. This semaphore caps how many HTTP requests can be inside
# the pipeline at once, so concurrent requests don't oversubscribe shared
# cores by each spawning their own cpu_count()-sized worker pool.
#
# Created lazily (not at import time) because asyncio.Semaphore needs a
# running event loop in some Python versions; a FastAPI dependency function
# creates and caches it on first use instead.
_semaphore: Optional[asyncio.Semaphore] = None


def get_semaphore() -> asyncio.Semaphore:
    global _semaphore
    if _semaphore is None:
        _semaphore = asyncio.Semaphore(get_settings().CONCURRENCY_LIMIT)
    return _semaphore
