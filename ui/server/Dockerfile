FROM hera-base:latest

COPY server/requirements.txt .

RUN pip install --no-cache-dir -r requirements.txt

COPY server/ .

EXPOSE 8000
EXPOSE 5678

CMD ["python", "-m", "debugpy", "--listen", "0.0.0.0:5678", "-m", "uvicorn", "main:app", "--host", "0.0.0.0", "--port", "8000"]
