all: ./src/sckmerdb_inspect.cpp ./src/sckmerdb_build.cpp ./src/gt_pro.cpp
	g++ -std=c++11 ./src/sckmerdb_inspect.cpp -o ./sckmerdb_inspect -O3
	g++ -std=c++11 ./src/sckmerdb_build.cpp -o ./sckmerdb_build -O3
	g++ -std=c++11 ./src/gt_pro.cpp -o ./gt_pro -O3 -pthread

clean:
	rm ./sckmerdb_inspect ./sckmerdb_build ./gt_pro

