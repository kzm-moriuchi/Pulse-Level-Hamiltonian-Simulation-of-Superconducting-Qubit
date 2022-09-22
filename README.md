# Pulse-Level-Hamiltonian-Simulation-of-Superconducting-Qubit
QuTiPのmaster equation solverである"mesolve"を用いて、T1,T2ノイズで特徴づけられる超伝導量子ビット系の量子操作をパルスレベルのハミルトニアンでシミュレートするプログラムを格納しています。

# Requirement
QuTiP 4.7.0 以上

# Installation
pip install qutip

# Usage
・main.ipynbの使い方  
上から順に実行すればOK。関数に与えるべきパラメータは関数内に記述しているので、パラメータを調整して好きなパルスをQubitに照射することができる。

・lib.H_gen_module.py  
今は、超伝導量子ビットの結合系を仮定しSchrieffer-Wolf変換後のハミルトニアンを生成するモジュールとしている。
