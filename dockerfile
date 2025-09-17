# Use Ubuntu 22.04 as the base image
FROM ubuntu:22.04

# Set environment variables to avoid interactive prompts during package installation
ENV DEBIAN_FRONTEND=noninteractive

# Install required packages for pyenv and building Python
RUN apt-get update && \
    apt-get install -y \
    curl \
    git \
    build-essential \
    libreadline-dev \
    libssl-dev \
    libbz2-dev \
    libsqlite3-dev \
    libffi-dev \
    zlib1g-dev && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Install pyenv
RUN curl https://pyenv.run | bash

# Set up environment variables for pyenv
ENV PATH="/root/.pyenv/bin:/root/.pyenv/shims:${PATH}"

# Install Python 3.9.13 using pyenv
RUN pyenv install 3.9.13 && \
    pyenv global 3.9.13

# Set the working directory
WORKDIR /app

# # Copy the requirements file into the container
# COPY requirements.txt .

# # Install pip for the installed Python version
# RUN python -m ensurepip && \
#     python -m pip install --no-cache-dir --upgrade pip

# # Install the Python dependencies
# RUN python -m pip install --no-cache-dir -r requirements.txt

# # Copy the rest of the application code
# COPY . .

# # Command to run the application
# CMD ["python", "app.py"]
