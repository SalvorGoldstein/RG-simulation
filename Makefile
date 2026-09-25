ccc:
	g++ -std=c++17 -O2 \
  -I autodiff \
  -I eigen-3.4.0 \
  main.cpp \
  -o main && ./main && rm ./main && python screen.py

