# bm-simulations

<p align="center">
  <img width="460" height="300" src="https://upload.wikimedia.org/wikipedia/commons/c/c2/Brownian_motion_large.gif">
</p>

## Simulations
- **Standard** Brownian Motion
- **Arithmetic** Brownian Motion (with drift and absorbing barrier)
- Geometric Brownian Motion (with ***absorbing barrier***)

## Running Docker Container
To keep simulations reproducible it is recommended they be run inside a Docker container. 

### Install Docker

See [Docker](https://docs.docker.com/get-docker/) for instructions. 

### Build Image

In **terminal**, change into the directory for this repo and run the following command: 

    docker build .

### Run Container
Once the container has been built, check the **IMAGE ID** by running: 

    docker images

Using your latest image ID run the following: 

    docker run -p 8888:8888 <your image id>

### Open Jupyter Notebook
Open your browser to the **URL** provided from the previous run command. It should look something like this: 
    
    http://127.0.0.1:8888/?token=26d71f7ec2d5c2f2deadaee53c36b8c7ef9a5918735b555f

This opens a Jupyter notebook. You should see the following: 

![enter image description here](http://collaboratescience.com/stack/random/jup.png)
Then open a new R notebook: 

![enter image description here](https://collaboratescience.com/stack/random/new_r.png)

You can now source all functions available inside bm_simulations.R

    source("bm-simulations/bm_simulations.R")
