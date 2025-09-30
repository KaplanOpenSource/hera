FROM hera:latest

WORKDIR /app/ui/server

COPY requirements.txt .

RUN pip install --no-cache-dir -r requirements.txt

COPY . .

EXPOSE 8000
EXPOSE 5678

# initializing with default project:
RUN sh ../../hera/scripts/run_mongo.sh && sh ../../hera/scripts/add_repo.sh

# TODO: run_mongo should be on bashrc, find out why it is not running.
CMD ["bash", "-c", "sh ../../hera/scripts/run_mongo.sh && python main.py"]
# CMD ["python", "-m", "debugpy", "--listen", "0.0.0.0:5678", "-m", "uvicorn", "main:app", "--host", "0.0.0.0", "--port", "8000"]
