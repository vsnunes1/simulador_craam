# simulador_craam

docker build --tag simulador_craam .
docker run -p 80:80 --name simulador -d simulador_craam