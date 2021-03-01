all: gt_pro mk_pool db_val db_dump db_uniq db_build
	@echo "GT-Pro build completed."

gt_pro:  ./src/gt_pro.cpp Makefile
	g++ -std=c++11 ./src/gt_pro.cpp -o ./bin/gt_pro -O3 -pthread

mk_pool: ./src/mk_pool.cpp Makefile
	g++ -std=c++11 ./src/mk_pool.cpp -o ./bin/mk_pool -O3

db_val: ./src/db_val.cpp Makefile
	g++ -std=c++11 ./src/db_val.cpp -o ./bin/db_val -O3 -pthread

db_dump: ./src/db_dump.cpp Makefile
	g++ -std=c++11 ./src/db_dump.cpp -o ./bin/db_dump -O3

db_uniq: ./src/db_uniq.cpp Makefile
	g++ -std=c++11 ./src/db_uniq.cpp -o ./bin/db_uniq -O3 

db_build: src/db_build.cpp Makefile
	g++ -std=c++11 ./src/db_build.cpp -o ./bin/db_build -O3 -pthread

clean:
	rm ./bin/db_build ./bin/gt_pro ./bin/mk_pool ./bin/db_val ./bin/db_dump ./bin/db_uniq

reformat:
	clang-format -style="{BasedOnStyle: llvm, ColumnLimit: 128}" src/sckmerdb_build.cpp > tmp-1.cpp && mv tmp-1.cpp src/sckmerdb_build.cpp
	clang-format -style="{BasedOnStyle: llvm, ColumnLimit: 128}" src/gt_pro.cpp > tmp-2.cpp && mv tmp-2.cpp src/gt_pro.cpp
