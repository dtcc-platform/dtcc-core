import os, pathlib, sys
import asyncio
import logging
import uvicorn

project_dir = str(pathlib.Path(__file__).resolve().parents[1])
sys.path.append(project_dir)

from src.app_fast import app as app_fastapi
from src.task_scheduler_publisher import scheduler as app_rocketry

fastapi_port = int(os.environ.get("FASTAPI_PORT", "8090"))

class Server(uvicorn.Server):
    """Customized uvicorn.Server
    
    Uvicorn server overrides signals and we need to include
    Rocketry to the signals."""
    def handle_exit(self, sig: int, frame) -> None:
        app_rocketry.session.shut_down()
        return super().handle_exit(sig, frame)


async def main():
    "Run Rocketry and FastAPI"
    server = Server(config=uvicorn.Config(app_fastapi, workers=1, loop="asyncio", port=fastapi_port, host="0.0.0.0",log_level='info'))

    api = asyncio.create_task(server.serve())
    sched = asyncio.create_task(app_rocketry.serve())

    await asyncio.wait([sched, api])

if __name__ == "__main__":
    # Print Rocketry's logs to terminal
    logger = logging.getLogger("rocketry.task")
    logger.addHandler(logging.StreamHandler())

    # Run both applications
    asyncio.run(main())