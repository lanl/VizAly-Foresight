# Travis CI and Docker for VizAly-CBench

## Build Docker image

To build docker image and tag as ``v0``, set ``${DOCKERHUB_USERNAME}`` to your DockerHub username and do
```
TAG=v0
docker build -t cbench .
docker tag cbench:latest ${DOCKERHUB_USERNAME}/vizaly-cbench:${TAG}
docker push ${DOCKERHUB_USERNAME}/vizaly-cbench:${TAG}
```

To use the new image in Travis CI update ``.travis.yml`` at the top-level of this repository to use your new tagged Docker image.
