ccc: main.cpp screen.py
	g++ -std=c++17 main.cpp -I autodiff -I eigen-3.4.0 -o main && ./main &&rm ./main && python screen.py
 

