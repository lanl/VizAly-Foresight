# Travis CI and Docker for VizAly-CBench

## Build Docker image

To build docker image and tag as ``v0``, if needed set your proxy (``${http_proxy}`` and ``${https_proxy}``), set ``${DOCKERHUB_USERNAME}`` to your DockerHub username, and then do:
```
TAG=v0
docker build -t cbench --build-arg http_proxy=${http_proxy} --build-arg https_proxy=${https_proxy} .
docker tag cbench:latest ${DOCKERHUB_USERNAME}/vizaly-cbench:${TAG}
docker push ${DOCKERHUB_USERNAME}/vizaly-cbench:${TAG}
```

To use the new image in Travis CI update ``.travis.yml`` at the top-level of this repository to use your new tagged Docker image.
