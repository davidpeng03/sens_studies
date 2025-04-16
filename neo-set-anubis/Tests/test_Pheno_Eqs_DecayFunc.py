import unittest
import sys, os

sys.path.append(os.getcwd())
import numpy as np
import scipy.integrate as integrate
from Core.QCDRunner import QCDRunner
from Core.Paramaters import SimulationParameters
from Pipeline.HNL.Equations.Pheno_Eqs import DecayFunctions, xSqrt, lambda_func, IFunc, Int_func, Lfunc

class TestDecayFunctions(unittest.TestCase):

    def setUp(self):
        self.MX = 100.0
        self.Ve = 0.1
        self.Vmu = 0.1
        self.Vtau = 0.1
        self.decay_functions = DecayFunctions(self.MX, self.Ve, self.Vmu, self.Vtau)

    def test_xSqrt(self):
        self.assertAlmostEqual(xSqrt(0.5), np.sqrt(1 - 4 * 0.5**2))

    def test_lambda_func(self):
        self.assertEqual(lambda_func(1, 2, 3), 1**2 + 2**2 + 3**2 - 2 * (1*2 + 1*3 + 2*3))

    def test_IFunc(self):
        result = IFunc(0.5, 0.3, 0.2, 0.1)
        expected = (1/0.5) * (0.5 - 0.1**2 - 0.2**2) * (1 + 0.3**2 - 0.5) * np.sqrt(lambda_func(0.5, 0.1**2, 0.2**2) * lambda_func(1, 0.5, 0.3**2))
        self.assertAlmostEqual(result, expected)

    def test_Int_func(self):
        result = Int_func(0.3, 0.2, 0.1)
        expected = 12 * integrate.quad(IFunc, (0.2 + 0.1)**2, (1 - 0.3)**2, args=(0.3, 0.2, 0.1))[0]
        self.assertAlmostEqual(result, expected)

    def test_Lfunc(self):
        self.assertAlmostEqual(Lfunc(0.5), np.log((1 - 3*0.5**2 - (1 - 0.5**2) * xSqrt(0.5)) / (0.5**2 * (1 + xSqrt(0.5)))))

    def test_is_approximately_equal(self):
        self.assertTrue(self.decay_functions.is_approximately_equal(1.0, 1.01, 0.01))
        self.assertFalse(self.decay_functions.is_approximately_equal(1.0, 1.02, 0.01))

    def test_set_get_MX(self):
        self.decay_functions.set_MX(200.0)
        self.assertEqual(self.decay_functions.get_MX(), 200.0)

    def test_set_get_couplings(self):
        self.decay_functions.set_couplings(0.2, 0.3, 0.4)
        self.assertEqual(self.decay_functions.get_couplings(), {11: 0.2, 13: 0.3, 15: 0.4})

    def test_dec_lep(self):
        result = self.decay_functions.dec_lep()
        self.assertEqual(result, 42)  # Placeholder value, replace with actual expected result

    def test_dec_invis(self):
        result = self.decay_functions.dec_invis()
        self.assertEqual(result, 42)  # Placeholder value, replace with actual expected result

    def test_dec_meson_v(self):
        result = self.decay_functions.dec_meson_v()
        self.assertEqual(result, 42)  # Placeholder value, replace with actual expected result

    def test_dec_meson_lep(self):
        result = self.decay_functions.dec_meson_lep()
        self.assertEqual(result, 42)  # Placeholder value, replace with actual expected result

    def test_QCDcorr(self):
        result = self.decay_functions.QCDcorr()
        self.assertEqual(result, 42)  # Placeholder value, replace with actual expected result

    def test_vTotHad(self):
        result = self.decay_functions.vTotHad()
        self.assertEqual(result, 42)  # Placeholder value, replace with actual expected result

    def test_vHad(self):
        result = self.decay_functions.vHad()
        self.assertEqual(result, 42)  # Placeholder value, replace with actual expected result

    def test_lHad(self):
        result = self.decay_functions.lHad()
        self.assertEqual(result, 42)  # Placeholder value, replace with actual expected result

    def test_lTotHad(self):
        result = self.decay_functions.lTotHad()
        self.assertEqual(result, 42)  # Placeholder value, replace with actual expected result

    def test_lAUD(self):
        result = self.decay_functions.lAUD(1.0, 0.1, 1.0, 1.0, 0.1)
        self.assertEqual(result, 42)  # Placeholder value, replace with actual expected result

    def test_vAff(self):
        result = self.decay_functions.vAff(1.0, 0.1)
        self.assertEqual(result, 42)  # Placeholder value, replace with actual expected result

    def test_vff(self):
        result = self.decay_functions.vff(1.0)
        self.assertEqual(result, 42)  # Placeholder value, replace with actual expected result

    def test_lHad_indi(self):
        result = self.decay_functions.lHad_indi(1.0, 0.1)
        self.assertEqual(result, 42)  # Placeholder value, replace with actual expected result

    def test_Nus(self):
        result = self.decay_functions.Nus()
        self.assertEqual(result, 42)  # Placeholder value, replace with actual expected result

    def test_lAhP(self):
        result = self.decay_functions.lAhP(1.0, 1.0, 1.0, 0.1, 0.1)
        self.assertEqual(result, 42)  # Placeholder value, replace with actual expected result

    def test_vAhP(self):
        result = self.decay_functions.vAhP(1.0, 1.0)
        self.assertEqual(result, 42)  # Placeholder value, replace with actual expected result

    def test_vhP(self):
        result = self.decay_functions.vhP(1.0, 1.0)
        self.assertEqual(result, 42)  # Placeholder value, replace with actual expected result

    def test_lAhV(self):
        result = self.decay_functions.lAhV(1.0, 1.0, 1.0, 0.1, 0.1)
        self.assertEqual(result, 42)  # Placeholder value, replace with actual expected result

    def test_vAhV(self):
        result = self.decay_functions.vAhV(1.0, 1.0, 1.0)
        self.assertEqual(result, 42)  # Placeholder value, replace with actual expected result

    def test_vhV(self):
        result = self.decay_functions.vhV(1.0, 1.0, 1.0, 0.1, 0.1, 0.1)
        self.assertEqual(result, 42)  # Placeholder value, replace with actual expected result

    def test_dec_tot(self):
        result = self.decay_functions.dec_tot()
        self.assertEqual(result, 42)  # Placeholder value, replace with actual expected result

if __name__ == "__main__":
    unittest.main()
