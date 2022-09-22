"""ハミルトニアンの係数(関数)を生成するモジュール"""

#モジュールのインポート
import numpy as np
from qutip import tensor, sigmap, sigmam, sigmax, sigmay, sigmaz, qeye, mesolve, sesolve, basis, identity

#クラス定義
class H_gen:

    #コンストラクタ
    def __init__(self, nqubits, pair, g, d_p, d_m, w, w_tilde, xi, nu_p, nu_m, mu_p, mu_m, phi, Omega):
        #初期生成変数
        self.nqubits = nqubits
        self.pair = pair
        self.g = g
        self.d_p = d_p
        self.d_m = d_m
        self.w = w
        self.w_tilde = w_tilde
        self.xi = xi
        self.nu_p = nu_p
        self.nu_m = nu_m
        self.mu_p = mu_p
        self.mu_m = mu_m
        self.phi = phi
        self.Omega = Omega

        #ハミルトニアンのパウリ表現を定義
        #qubit 0 → 1へのパルス
        #XIの項
        self.xi_0 = [qeye(2) for _ in range(self.nqubits)]
        self.xi_0[self.pair[0]] = sigmax()
        self.XI_0 = tensor(self.xi_0)
        #YIの項
        self.yi_0 = [qeye(2) for _ in range(self.nqubits)]
        self.yi_0[self.pair[0]] = sigmay()
        self.YI_0 = tensor(self.yi_0)
        #IXの項
        self.ix_0 = [qeye(2) for _ in range(self.nqubits)]
        self.ix_0[self.pair[1]] = sigmax()
        self.IX_0 = tensor(self.ix_0)
        #IYの項
        self.iy_0 = [qeye(2) for _ in range(self.nqubits)]
        self.iy_0[self.pair[1]] = sigmay()
        self.IY_0 = tensor(self.iy_0)
        #ZXの項(若い方のqubitからのドライブ)
        self.zx_0 = [qeye(2) for _ in range(self.nqubits)]
        self.zx_0[self.pair[0]], self.zx_0[self.pair[1]] = sigmaz(), sigmax()
        self.ZX_0 = tensor(self.zx_0)
        #ZYの項(若い方のqubitからのドライブ)
        self.zy_0 = [qeye(2) for _ in range(self.nqubits)]
        self.zy_0[self.pair[0]], self.zy_0[self.pair[1]] = sigmaz(), sigmay()
        self.ZY_0 = tensor(self.zy_0)

        #qubit 1 → 0へのパルス
        #XIの項
        self.xi_1 = [qeye(2) for _ in range(self.nqubits)]
        self.xi_1[self.pair[1]] = sigmax()
        self.XI_1 = tensor(self.xi_1)
        #YIの項
        self.yi_1 = [qeye(2) for _ in range(self.nqubits)]
        self.yi_1[self.pair[1]] = sigmay()
        self.YI_1 = tensor(self.yi_1)
        #IXの項
        self.ix_1 = [qeye(2) for _ in range(self.nqubits)]
        self.ix_1[self.pair[0]] = sigmax()
        self.IX_1 = tensor(self.ix_1)
        #IYの項
        self.iy_1 = [qeye(2) for _ in range(self.nqubits)]
        self.iy_1[self.pair[0]] = sigmay()
        self.IY_1 = tensor(self.iy_1)
        #ZXの項(若い方のqubitからのドライブ)
        self.zx_1 = [qeye(2) for _ in range(self.nqubits)]
        self.zx_1[self.pair[1]], self.zx_1[self.pair[0]] = sigmaz(), sigmax()
        self.ZX_1 = tensor(self.zx_1)
        #ZYの項(若い方のqubitからのドライブ)
        self.zy_1 = [qeye(2) for _ in range(self.nqubits)]
        self.zy_1[self.pair[1]], self.zy_1[self.pair[0]] = sigmaz(), sigmay()
        self.ZY_1 = tensor(self.zy_1)


    #ハミルトニアンの係数(関数)を定義
    #qubit 0 → 1へのパルス
    def XI_coef_0(self, t, args):
        return 1/2 * self.Omega[self.pair[0]] * np.cos((self.w[self.pair[0]] - self.w_tilde[self.pair[0]]) * t + self.phi[self.pair[0]])

    def YI_coef_0(self, t, args):
        return -1/2 * self.Omega[self.pair[0]] * np.sin((self.w[self.pair[0]] - self.w_tilde[self.pair[0]]) * t + self.phi[self.pair[0]])

    def IX_coef_0(self, t, args):
        return -self.nu_p[self.pair[0]]/2 * self.Omega[self.pair[0]] * np.cos((self.w[self.pair[0]] - self.w_tilde[self.pair[1]]) * t + self.phi[self.pair[0]])

    def IY_coef_0(self, t, args):
        return self.nu_p[self.pair[0]]/2 * self.Omega[self.pair[0]] * np.sin((self.w[self.pair[0]] - self.w_tilde[self.pair[1]]) * t + self.phi[self.pair[0]])

    def ZX_coef_0(self, t, args):
        return -self.mu_p[self.pair[0]]/2 * self.Omega[self.pair[0]] * np.cos((self.w[self.pair[0]] - self.w_tilde[self.pair[1]]) * t + self.phi[self.pair[0]])

    def ZY_coef_0(self, t, args):
        return self.mu_p[self.pair[0]]/2 * self.Omega[self.pair[0]] * np.sin((self.w[self.pair[0]] - self.w_tilde[self.pair[1]]) * t + self.phi[self.pair[0]])


    #qubit 1 → 0へのパルス
    def XI_coef_1(self, t, args):
        return 1/2 * self.Omega[self.pair[1]] * np.cos((self.w[self.pair[1]] - self.w_tilde[self.pair[1]]) * t + self.phi[self.pair[1]])

    def YI_coef_1(self, t, args):
        return -1/2 * self.Omega[self.pair[1]] * np.sin((self.w[self.pair[1]] - self.w_tilde[self.pair[1]]) * t + self.phi[self.pair[1]])

    def IX_coef_1(self, t, args):
        return -self.nu_m[self.pair[0]]/2 * self.Omega[self.pair[1]] * np.cos((self.w[self.pair[1]] - self.w_tilde[self.pair[0]]) * t + self.phi[self.pair[1]])

    def IY_coef_1(self, t, args):
        return self.nu_m[self.pair[0]]/2 * self.Omega[self.pair[1]] * np.sin((self.w[self.pair[1]] - self.w_tilde[self.pair[0]]) * t + self.phi[self.pair[1]])

    def ZX_coef_1(self, t, args):
        return -self.mu_m[self.pair[0]]/2 * self.Omega[self.pair[1]] * np.cos((self.w[self.pair[1]] - self.w_tilde[self.pair[0]]) * t + self.phi[self.pair[1]])

    def ZY_coef_1(self, t, args):
        return self.mu_m[self.pair[0]]/2 * self.Omega[self.pair[1]] * np.sin((self.w[self.pair[1]] - self.w_tilde[self.pair[0]]) * t + self.phi[self.pair[1]])


    #全ハミルトニアンを係数とパウリ演算子のリストで返す関数
    def H_td_list(self):
        # self.H_td = [[self.ZZ, self.ZZ_coef], [self.XI_0, self.XI_coef_0], [self.YI_0, self.YI_coef_0], [self.IX_0 ,self.IX_coef_0], [self.IY_0, self.IY_coef_0], [self.ZX_0, self.ZX_coef_0], [self.ZY_0, self.ZY_coef_0], [self.XI_1, self.XI_coef_1], [self.YI_1, self.YI_coef_1], [self.IX_1 ,self.IX_coef_1], [self.IY_1, self.IY_coef_1], [self.ZX_1, self.ZX_coef_1], [self.ZY_1, self.ZY_coef_1]]
        self.H_td = [[self.XI_0, self.XI_coef_0], [self.YI_0, self.YI_coef_0], [self.IX_0 ,self.IX_coef_0], [self.IY_0, self.IY_coef_0], [self.ZX_0, self.ZX_coef_0], [self.ZY_0, self.ZY_coef_0], [self.XI_1, self.XI_coef_1], [self.YI_1, self.YI_coef_1], [self.IX_1 ,self.IX_coef_1], [self.IY_1, self.IY_coef_1], [self.ZX_1, self.ZX_coef_1], [self.ZY_1, self.ZY_coef_1]]
        return self.H_td

    def check(self):
        return self.w[self.dep][self.pair[0]] - self.w_tilde[self.pair[0]]


