all: sckmerdb_inspect sckmerdb_build gtpro
	@echo "GTPro build completed."

gtpro:  ./src/gt_pro.cpp Makefile
	g++ -std=c++11 ./src/gt_pro.cpp -o ./gt_pro -O3 -pthread

sckmerdb_build: src/sckmerdb_build.cpp Makefile
	g++ -std=c++11 ./src/sckmerdb_build.cpp -o ./sckmerdb_build -O3 -pthread

clean:
	rm ./sckmerdb_build ./gt_pro
