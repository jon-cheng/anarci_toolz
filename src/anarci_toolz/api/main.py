import uvicorn
from fastapi import FastAPI

from .routes import router

app = FastAPI(title="anarci-toolz API")
app.include_router(router)


def run() -> None:
    uvicorn.run("anarci_toolz.api.main:app", host="0.0.0.0", port=8000)


if __name__ == "__main__":
    run()
