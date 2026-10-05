ccc: main.cpp Schwartzschield_aff.py
	g++ -std=c++17 main.cpp outfunc.cpp -I autodiff -I eigen-3.4.0 -o main && ./main &&rm ./main && python Schwartzschield_aff.py

