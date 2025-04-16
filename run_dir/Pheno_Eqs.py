import numpy as np
import scipy.integrate as integrate

from common_funcs import * 
from const_params import *

# ------- Functions
def IFunc(x,xu,xd,xl):
    return (1/x)*(x-xl**2-xd**2)*(1+xu**2-x)*np.sqrt(lambda_func(x,xl**2,xd**2)*lambda_func(1,x,xu**2))

def Int_func(xu,xd,xl):
    return 12*integrate.quad(IFunc, (xd + xl)**2, (1- xu)**2, args=(xu,xd,xl))[0]

# ------- Parameters
C1fList = [ (1/4)*(1 - (8/3)*sw2 + (32/9)*sw2**2), (1/4)*(1 - (4/3)*sw2 + (8/9)*sw2**2), (1/4)*(1 - 4*sw2 +8*sw2**2), (1/4)*(1 + 4*sw2 + 8*sw2**2) ]
C2fList = [ (1/3)*sw2*((4/3)*sw2 - 1), (1/6)*sw2*((2/3)*sw2 - 1), (1/2)*sw2*(2*sw2 - 1), (1/2)*sw2*(2*sw2 + 1) ]

Vlist = ["Ve","Vmu","Vtau"]

# ------- Pheno equations
# from 1805.08567

# np.heaviside(MX - Ml - Mh,1) make this into a function

# ---- Three Body

# Charged Current-Mediated Decays
# lA nuB lB & lA U D-bar (neutrino and fermion pair, ferimons different)
def lAUD(MX,Ml,Vl,MUp,MDown,Vij): 
    if MX >= (Ml + MUp + MDown):
        P1 = (Gf**2*MX**5)/(192*pi**3) 
        P2 = Int_func(MUp/MX, MDown/MX, Ml/MX)

        if (MUp in mQlist) and (MDown in mQlist): # quarks
            NW = 3*np.abs(Vij)**2 #Nc = 3 for quarks
        elif (MUp ==0 ) and (MDown in mLlist): # leptons
            NW = 1
        else: #otherwise
            NW = 0

        return (NW * P1 * np.abs(Vl)**2 * P2)
    else:
        return 0

# Neutral Current-Mediates Decays
# neutrino and fermion pair (same fermions)
def vAff(MX,Mf,Vi):
    if MX >= 2*Mf:
        xf = Mf/MX

        P1 = (Gf**2*MX**5)/(192*pi**3)

        ix = 10 # will cause error for wrong mass
        if (Mf in mUlist):
            ix = 0
        elif (Mf in mDlist):
            ix = 1
        elif (Mf in mLlist):
            if ( Vi == 1 ): # same neutrino and lepton flavour
                ix = 3
            else:
                ix = 2

        if ix == 10:
            return 0
        else:
            C1f = C1fList[ix]
            C2f = C2fList[ix]

        P2 = (
            C1f * ( (1- 14*xf**2 - 2*xf**4 - 12*xf**6)*np.sqrt(1-4*xf**2)
                + 12*xf**4 * (xf**4 - 1) *Lfunc(xf))
            + 4*C2f *(xf**2 * (2+ 10*xf**2 - 12*xf**4)*xSqrt(xf)
                    + 6*xf**4 * (1 - 2*xf**2 + 2*xf**4)*Lfunc(xf))       
            )

        if (Mf in mQlist): # quarks
            NZ = 3 #Nc = 3 for quarks
        elif (Mf in mLlist): # leptons
            NZ = 1
        else: #otherwise
            NZ = 0

        return NZ * P1 * P2
    else:
        return 0

def vff(MX,Mf,Ve,Vmu,Vtau):
    outRes = 0
    VvalList = [Ve,Vmu,Vtau]

    for VA in Vlist:
        VAval = VvalList[Vlist.index(VA)]
        Vi = 0
        # if lepton, difference to whether neutrino match the lepton flavour
        if (Mf in mLlist):
            ix = np.where(np.round(mLlist,6) == np.round(Mf,6))[0][0]
            VB = Vlist[ix]

            if VA == VB:
                Vi = 1

        outRes += np.abs(VAval)**2 *vAff(MX,Mf,Vi)

    return outRes


# Three neutrinos
def Nus(MX,Ve,Vmu,Vtau):
    # (1 + dirac_delta(A,B)) summed over 3, where one is 2 => 1 + 1 + 2 = 4
    P1 = (Gf**2)/(192*pi**3)
    P2 = MX**5

    return P1 * (np.abs(Ve)**2 + np.abs(Vmu)**2 + np.abs(Vtau)**2) * P2

# ---- Two Body
# Charged Pseudo-Scalar
def lAhP(MX, Ml, Mh, fh, Vud, VA):
    if MX >= Ml + Mh:
        xl = Ml/MX
        xh = Mh/MX

        P1 = (Gf**2 *fh**2 * np.abs(Vud)**2)/(16*pi)
        P2 = MX**3 * ( (1-xl**2)**2 - xh**2*(1+xl**2) )*np.sqrt(lambda_func(1,xh**2,xl**2))

        return np.abs(VA)**2 * P1 * P2
    else:
        return 0

# Neutral Pseudo-Scalar
def vAhP(MX, Mh, fh):
    xh = Mh/MX

    P1 = (Gf**2 *fh**2)/(32*pi)
    P2 = MX**3 * (1- xh**2)**2

    return P1 * P2

def vhP(MX, Mh, fh, Ve, Vmu, Vtau):
    if MX >= Mh:
        out_val = (np.abs(Ve)**2 * vAhP(MX, Mh, fh) 
            + np.abs(Vmu)**2 * vAhP(MX, Mh, fh) 
            + np.abs(Vtau)**2 * vAhP(MX, Mh, fh))
        return out_val
    else:
        return 0

# Charged Vector
def lAhV(MX, Ml, Mh, gh, Vud, VA):
    if MX >= Ml + Mh:
        xl = Ml/MX
        xh = Mh/MX

        P1 = (Gf**2 *gh**2 * np.abs(Vud)**2)/(16*pi * Mh**2)
        P2 = MX**3 * ( (1-xl**2)**2 + xh**2*(1+xl**2) -2*xh**4 )*np.sqrt(lambda_func(1,xh**2,xl**2))

        return np.abs(VA)**2 * P1 * P2
    else:
        return 0

# Neutral Vector
def vAhV(MX, Mh, kh, gh):
    xh = Mh/MX

    P1 = (Gf**2 * kh**2 *gh**2)/(32*pi + Mh**2)
    P2 = MX**3 * (1 + 2*xh**2) * (1- xh**2)**2

    return P1 * P2

def vhV(MX, Mh, kh, gh,Ve,Vmu,Vtau):
    if MX >= Mh:
        out_val = (np.abs(Ve)**2 * vAhV(MX, Mh, kh, gh) 
            + np.abs(Vmu)**2 * vAhV(MX, Mh, kh, gh) 
            + np.abs(Vtau)**2 * vAhV(MX, Mh, kh, gh))
        return out_val
    else:
        return 0

# ---- Total Hadronic
# MX > 2*mPion

def QCDcorr(MX):
    # options: aSrun(MX), aSrun_1(MX), alpha_vals[ix]
    ix = np.where(np.round(MX,1) == m_vals)[0][0]
    aS_val = alpha_vals[ix]
    
    return (aS_val/pi) + 5.2*(aS_val/pi)**2 + 26.4*(aS_val/pi)**3

def vTotHad(MX,Ve,Vmu,Vtau):
    dec_out = 0

    for i in range(0,len(mQlist)-1):
        Mq = mQlist[i]
        M_test = ["mKaon","mD","mB",]

        if i <= 1:
            dec_out += vff(MX,Mq,Ve,Vmu,Vtau)
        else:
            mMes = mMasses[M_test[i-2]]
            if MX >= 2*mMes:
                dec_out += np.sqrt(1-4*(mMes**2)/(MX**2))*vff(MX,Mq,Ve,Vmu,Vtau)

    return (1+QCDcorr(MX))*dec_out

def lTotHad(MX,Ve,Vmu,Vtau):
    dec_out = 0

    for i in range(0,len(mLlist)):
        Ml = mLlist[i]
        """
        match i:
                case 0:
                    Vl = Ve
                case 1:
                    Vl = Vmu
                case 2:
                    Vl = Vtau
                    """
               # Replace match-case with if-elif
        if i == 0:
            Vl = Ve
        elif i == 1:
            Vl = Vmu
        elif i == 2:
            Vl = Vtau

        for j in range(0,len(mUlist)):
            for k in range(0,len(mDlist)):
                Vij = CKM_matrix[j][k]

                MUp = mUlist[j]
                MDown = mUlist[k]

                dec_out += lAUD(MX,Ml,Vl,MUp,MDown,Vij)

        #dec_out += lAUD(MX,Ml,Vl,mUlist[0],mDlist[0],CKM_matrix[0][0])
        #dec_out += lAUD(MX,Ml,Vl,mUlist[0],mDlist[1],CKM_matrix[0][1])

    return (1+QCDcorr(MX))*dec_out

# ------- Decay Widths ---------
# ------- Mesons
# Sum over final states with neutrino
def dec_meson_v(MX, Ve, Vmu, Vtau):
    dec_out = 0

    # Neutral Pseudo-Scalars
    for i in range(0,len(PseMes0List)):
        Mh = mMasses["m"+PseMes0List[i]] # mass
        fh = mDecays["f"+PseMes0List[i]] # decay constant

        dec_out += vhP(MX,Mh,fh,Ve,Vmu,Vtau)

    # Neutral Vectors
    for i in range(0,len(VecMes0List)):
        Mh = mMasses["m"+VecMes0List[i]] # mass
        kh = mKh["k"+VecMes0List[i]] # constant
        gh = mgh["g"+VecMes0List[i]]

        dec_out += vhV(MX, Mh, kh, gh,Ve,Vmu,Vtau)

    return dec_out

# Sum over lepton final states
def dec_meson_lep(MX, Ve, Vmu, Vtau):
    dec_out = 0

    for i in range(0,len(PseMesList)):
        Mh = mMasses["m"+PseMesList[i]] # mass
        fh = mDecays["f"+PseMesList[i]] # decay constant
        Vud = mCKM["V"+PseMesList[i]] # CKM matrix

        for j in range(0,len(mLlist)):
            Ml = mLlist[j]
            # match j:
            #     case 0:
            #         VA = Ve
            #     case 1:
            #         VA = Vmu
            #     case 2:
            #         VA = Vtau

            if j == 0:
                VA = Ve
            elif j == 1:
                VA = Vmu
            elif j == 2:
                VA = Vtau

            dec_out += lAhP(MX, Ml, Mh, fh, Vud, VA)

    for i in range(0,len(VecMesList)):
        Mh = mMasses["m"+VecMesList[i]] # mass
        gh = mgh["g"+VecMes0List[i]] # decay constant
        Vud = mCKM["V"+VecMesList[i]] # CKM matrix

        for j in range(0,len(mLlist)):
            Ml = mLlist[j]

            # match j:
            #     case 0:
            #         VA = Ve
            #     case 1:
            #         VA = Vmu
            #     case 2:
            #         VA = Vtau

            if j == 0:
                VA = Ve
            elif j == 1:
                VA = Vmu
            elif j == 2:
                VA = Vtau
            
            dec_out += lAhV(MX, Ml, Mh, gh, Vud, VA)
    return dec_out

# ------- Leptonic
def dec_lep(MX, Ve, Vmu, Vtau):
    dec_out = 0

    for i in range(0,len(mLlist)):
        Ml = mLlist[i]
        
        # match i:
        #         case 0:
        #             Vl = Ve
        #         case 1:
        #             Vl = Vmu
        #         case 2:
        #             Vl = Vtau
        if i == 0:
            Vl = Ve
        elif i == 1:
            Vl = Vmu
        elif i == 2:
            Vl = Vtau

        for j in range(0,len(mLlist)):
            Vij = 0 # only for quarks

            MUp = 0 # for neutrino
            MDown = mLlist[j]

            dec_out += lAUD(MX,Ml,Vl,MUp,MDown,Vij)

        dec_out += vff(MX,Ml,Ve,Vmu,Vtau)

    return dec_out

# ------- Invisible
def dec_invis(MX,Ve,Vmu,Vtau):
    return Nus(MX,Ve,Vmu,Vtau)

# ---- Multi-Hadron Final States
# MX > 2*mPion

def vHad(MX,Ve,Vmu,Vtau):
    if MX > 1:#2*mMasses["mPion"]:
        return (vTotHad(MX,Ve,Vmu,Vtau) - dec_meson_v(MX, Ve, Vmu, Vtau))
    else:
        return 0
    
def lHad(MX,Ve,Vmu,Vtau):
    if MX > 1:#2*mMasses["mPion"]:
        return (lTotHad(MX,Ve,Vmu,Vtau) - dec_meson_lep(MX, Ve, Vmu, Vtau))
    else:
        return 0
    
def lHad_indi(MX,Ml,Vl):
    if MX > 1: #2*mMasses["mPion"]:
        dec_out = 0

        for j in range(0,len(mUlist)):
            for k in range(0,len(mDlist)):
                Vij = CKM_matrix[j][k]

                MUp = mUlist[j]
                MDown = mUlist[k]

                dec_out += lAUD(MX,Ml,Vl,MUp,MDown,Vij)

        for i in range(0,len(PseMesList)):
            Mh = mMasses["m"+PseMesList[i]]
            fh = mDecays["f"+PseMesList[i]]
            Vud = mCKM["V"+PseMesList[i]]

            dec_out -= lAhP(MX, Ml, Mh, fh, Vud, Vl)

        for i in range(0,len(VecMesList)):
            Mh = mMasses["m"+VecMesList[i]]
            gh = mgh["g"+VecMes0List[i]]
            Vud = mCKM["V"+VecMesList[i]]
                
            dec_out -= lAhV(MX, Ml, Mh, gh, Vud, Vl)

            return dec_out
    else:
        return 0

# ------- Total
def dec_tot(MX,Ve,Vmu,Vtau):
    return ( 
        dec_invis(MX,Ve,Vmu,Vtau)
        + dec_meson_v(MX,Ve,Vmu,Vtau)
        + dec_meson_lep(MX, Ve, Vmu, Vtau)
        + dec_lep(MX,Ve,Vmu,Vtau)
        + vHad(MX,Ve,Vmu,Vtau)
        + lHad(MX,Ve,Vmu,Vtau)
    )