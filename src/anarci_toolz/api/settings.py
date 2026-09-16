from functools import lru_cache

from pydantic_settings import BaseSettings


class Settings(BaseSettings):
    MAX_SEQUENCES: int = 20
    CONCURRENCY_LIMIT: int = 1
    SUBPROCESS_TIMEOUT_SECONDS: int = 60
    WORKER_COUNT: int = 1


@lru_cache
def get_settings() -> Settings:
    return Settings()
