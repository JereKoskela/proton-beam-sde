import math
import numpy as np
from endf_parserpy import EndfParser

#Script to output Elastic cross section data for hydrogen ready to be parsed into SDE model
#Cross sections are evaluated as described in https://iopscience.iop.org/article/10.1088/1361-6560/ae5586 DOI 10.1088/1361-6560/ae5586

#Cuts off cross section data at this forward scattering angle in RAD. Can be kept small as the variable parameter rutherford_cutoff in the main code
# can be used to modify this cutoff variably. 
cut_off = 0.0001
#Elastic scatter for hydrogen results in two protons produced. At large scatter angles the cross section data
#will generate the proton almost perpendicular to the original protons trajectory with almost zero energy
#back_scatter_cut_off_cm defines a cutoff for which the forward peaked proton is instead generated carrying
#almost all of the energy. Value of 0.9 corresponds to forward peaked proton carrying at least ~80% of total
#energy
back_scatter_cut_off_cm = 0.9
#Spherical brownian motion scattering density diverges from Moliere scattering density at 2.5x standard deviation
#Exclude cross section data above Gaussian_SD_cut_off's of this standard deviation. Should be at least 2.5x
#but can be fitted parameter
Gaussian_SD_cut_off = 2.5
## Timestep used for lynch_dahl small scattering SD, this should match the lynch_dahl_timestep
# in the config file
lynch_dahl_timestep = 0.05


# SD of multiple scattering contribution, individual scattering should start at 2.5SD. Taken with track length step 0.5mm (first constant in chi_c_sq).
def multiple_scattering_sd(e, z,a):
    mpcsq = 938.346  # mass of proton * speed of light squared, MeV
    pv = (2 * mpcsq + e) * e / (mpcsq + e)
    betasq = (2 * mpcsq + e) * e / ((mpcsq + e) ** 2)
    p = pv / math.sqrt(betasq) #momentum in MeV / c.
    chi_c_sq = lynch_dahl_timestep*z*(z+1)*0.157/(a*pv*pv)
    chi_a_sq_vec = 2.007e-5 * (z**(2 / 3)) * (1 + 3.34 * (z / (137 * math.sqrt(betasq)))**2) / (p * p)
    omega = chi_c_sq / chi_a_sq_vec
    F = 0.98
    v = omega / (2 * (1 - F))
    ret = math.sqrt(
        chi_c_sq * (((1 + v) * math.log(1 + v) / v) - 1) / (1 + F * F)
    )
    return ret

def Nuclear_trig_constants(n, eta, sign):  # recursively evaluate antiderivative
    if n == 0:
        return [0, 1]
    else:
        sin_const = (
            sign
            * (
                (n - 1) * Nuclear_trig_constants(n - 1, eta, sign)[0]
                + eta * Nuclear_trig_constants(n - 1, eta, sign)[1]
            )
            / ((n - 1) ** 2 + eta**2)
        )
        cos_const = (
            sign
            * (
                (n - 1) * Nuclear_trig_constants(n - 1, eta, sign)[1]
                - eta * Nuclear_trig_constants(n - 1, eta, sign)[0]
            )
            / ((n - 1) ** 2 + eta**2)
        )
        return [sin_const, cos_const]
#Test cases for Nuclear_trig_constants
#print("Testing Nuclear_trig_constants")
#for i in range(1,20):
#    test_input = i*random.uniform(0.0, 1.0)
#    print(abs(Nuclear_trig_constants(3, test_input, 1)[0]-(2-test_input**2)/((test_input**2+4)*(test_input**2+1)*test_input))<1e-10)
#    print(abs(Nuclear_trig_constants(3, test_input, 1)[1]+3/((test_input**2+4)*(test_input**2+1)))<1e-10)

def Trig_Integral_Eval(
    deg, cosval, sinval, y, eta, sign, dict
):  # evaluation of antiderivative at a given value y
    if deg in dict.keys():
        coeff = dict[deg]
    else:
        coeff = Nuclear_trig_constants(deg, eta, sign)
    return ((1 + (y * sign)) ** (deg - 1)) * (coeff[0] * sinval + coeff[1] * cosval)


def CDF_Poly_Eval(deg, eta, c, y, sign, dict):  # evaluation of integration by parts form
    cosval = np.cos(eta * np.log((1 + (y * sign)) / 2) + c)
    sinval = np.sin(eta * np.log((1 + (y * sign)) / 2) + c)
    temp = 0
    for i in range(deg + 1):
        temp += (
            Trig_Integral_Eval(i + 1, cosval, sinval, y, eta, sign, dict)
            * (y ** (deg - i))
            * math.factorial(deg)
            * ((-1) ** (i))
            / math.factorial(deg - i)
        )
    cosval0 = np.cos(eta * np.log(1 / 2) + c)
    sinval0 = np.sin(eta * np.log(1 / 2) + c)
    temp += (
        ((-1) ** (deg + 1))
        * math.factorial(deg)
        * Trig_Integral_Eval(deg + 1, cosval0, sinval0, 0, eta, sign, dict)
    )
    return temp


def Legendre_Coeff_Cal(deg):  # outputs coefficients of legendre polynomial as a list
    if deg == 0:
        return [1]
    elif deg == 1:
        return [0, 1]
    else:
        temp1 = Legendre_Coeff_Cal(deg - 1)
        temp1.insert(0, 0)
        for i in range(len(temp1)):
            temp1[i] = (2 * (deg - 1) + 1) * temp1[i] / deg
        temp2 = Legendre_Coeff_Cal(deg - 2)
        temp2.append(0)
        temp2.append(0)
        for i in range(len(temp2)):
            temp2[i] = (1 - deg) * temp2[i] / deg
        temp3 = []
        for i in range(len(temp2)):
            temp3.append(temp1[i] + temp2[i])
        return temp3
#Test cases for Legendre_Coeff_Cal
#print("Testing Legendre_Coeff_Cal")
#print(Legendre_Coeff_Cal(5) == [0,15/8,0,-70/8,0,63/8])
#print(Legendre_Coeff_Cal(8) == [35/128, 0, -1260/128, 0, 6930/128, 0, -12012/128, 0, 6435/128])

def Legendre_dict_production(deg):
    Legendre_dict = {}
    for i in range(deg):
        Legendre_dict[i] = Legendre_Coeff_Cal(i)
    return Legendre_dict
def Trig_dict_production(deg, eta, sign):
    Trig_dict = {}
    for i in range(deg):
        Trig_dict[i] = Nuclear_trig_constants(i, eta, sign)
    return Trig_dict
#Test cases for dict generators
#print("Testing Legendre_dict_production and Trig_dict_production")
#test_dict_Legendre = Legendre_dict_production(10)
#for key in test_dict_Legendre.keys():
#    print(test_dict_Legendre[key]==Legendre_Coeff_Cal(key))
#test_dict_Trig = Trig_dict_production(10, 1, 1)
#for key in test_dict_Trig.keys():
#    print(test_dict_Trig[key]==Nuclear_trig_constants(key, 1, 1))

def CDF_Legendre_Eval(
    deg, eta, c, y, sign, dict_leg, dict_trig
):  # evaluation of integration by parts for Legendre polynomial
    if deg in dict_leg.keys():
        Legendre_temp = dict_leg[deg]
    else:
        Legendre_temp = Legendre_Coeff_Cal(deg)
    out = 0
    for i in range(deg):
        out += Legendre_temp[i] * CDF_Poly_Eval(i, eta, c, y, sign, dict_trig)
    return out


def Legendre_int(deg, y, dict):  # evaluation of Legendre polynomial integral
    if deg in dict.keys():
        Legendre_temp = dict[deg]
    else:
        Legendre_temp = Legendre_Coeff_Cal(deg)
    Total = 0
    for i in range(len(Legendre_temp)):
        Total += Legendre_temp[i]*(y ** (i + 1)) / (i + 1)
    return Total


def Complex_Polar(re, im):
    r = np.sqrt(re**2 + im**2)
    theta = np.arctan2(im, re)
    return [r, theta]
#Evaluates eta as in Section 6.2.6 of the ENDF data manual
def eta_eval(Z, z, E, m):
    sqrtterm = 0.0496 * m / (E * 2)
    out = Z * z * np.sqrt(sqrtterm)
    return out

#Calculating Formula in A2 of https://iopscience.iop.org/article/10.1088/1361-6560/ae5586 DOI 10.1088/1361-6560/ae5586
def Total_Rate_calc(
    Z, z, E, A, a, y, two_s, a_coef_re, a_coef_im, b_coef, dict_trig_pv, dict_trig_nv, dict_leg
):  # total rate between 0 and y
    A_ratio = A / a
    Rutherford_Coef = ((1e-2) * 4.16 * ((Z * z) * (1 + A_ratio)) ** 2) / (
        (2 * A_ratio * E) ** 2
    )
    Rate = 0
    #Corrected error, sin eval at eta_val instead of 2*eta_val
    Rate += Rutherford_Coef * (
        (y / (1 - y**2))
        + ((((-1) ** two_s) / (two_s + 1)))
        * math.sin(eta_val * math.log((1 + y) / (1 - y)))
        / (2 * eta_val)
    )
    Nuclear_term_1 = 0
    for i in range(len(a_coef_im)):
        polar_temp = Complex_Polar(a_coef_re[i], a_coef_im[i])
        int_temp = 0
        int_temp += ((-1) ** i) * CDF_Legendre_Eval(i, eta_val, polar_temp[1], y, 1, dict_leg, dict_trig_pv)
        int_temp += CDF_Legendre_Eval(i, eta_val, polar_temp[1], y, -1, dict_leg, dict_trig_nv)
        Nuclear_term_1 += polar_temp[0] * (2 * i + 1) * int_temp / 2
    Nuclear_term_1 *= 2 * eta_val
    Rate -= Nuclear_term_1
    Nuclear_term_2 = 0
    for i in range(len(a_coef_im)):
        Nuclear_term_2 += (4 * i + 1) * b_coef[i] * Legendre_int(2 * i, y, dict_leg) / 2
    Rate += Nuclear_term_2
    Rate *= 2 * np.pi
    return Rate


def ErrorCheck(Data, sign):  # checks for monotonicity in CDF and CM to Lab frame
    if sign == 1:
        for i in range(1, len(Data)):
            if Data[i] < Data[i - 1]:
                print("Error in CDF monotonicity")
    if sign == -1:
        for i in range(1, len(Data)):
            if Data[i] >= Data[i - 1]:
                print("Error in CM to Lab monotonicity")
                print(Data[i], Data[i - 1])
    if sign != -1 and sign != 1:
        print("sign value error")
    return 0


def CM_to_Lab_Frame(ang0, E, A, a):
    ang = math.acos(ang0)
    mp = a * 938.346  # mass of proton * c^2, MeV
    mn = A * 938.346  # mass of colliding nucleus * c^2, MeV
    p = math.sqrt((E) * (E + 2 * mp))
    u = p / (E + mp + mn)  # typo should be 2mp
    g = 1 / math.sqrt(1 - u * u)
    e = E + mp
    v_ratio = u * (e - u * p) / (p - u * e)
    if np.fabs((g * (math.cos(ang) + v_ratio))) == 0:
        out = math.pi / 2
    else:
        out = math.atan(math.sin(ang) / (g * (math.cos(ang) + v_ratio)))
    if out < 0:
        out += math.pi
    return out


parser = EndfParser()
endf_dict = parser.parsefile("Splines/ENDF_data/hydrogen_el.txt")
density = {}
angle_discretize = np.linspace(0, 0.9, 100)
angle_discretize_end = []
for i in range(39):  # Setup interpolation points to match growth of rutherford crossection near boundary
    angle_discretize_end.append(np.sqrt((11 + 10 * i - 1) / (11 + 10 * i)))
for i in range(9000): #Sparser additional interpolation points if needed to minimise file size
    angle_discretize_end.append(np.sqrt((401 + 400 * i - 1) / (401 + 400 * i)))
angle_discretize_pos = np.concatenate((angle_discretize, angle_discretize_end))
Energy_discrete = []
#list of incident energy values
for i in endf_dict[6][2]["subsection"][1]["E"].keys():
    yy = []
    Total_rate = []
    Real_coef = []
    Imaginary_coef = []
    b_coef = []
    #incident energy value
    Energy = endf_dict[6][2]["subsection"][1]["E"][i] * 1e-6
    if Energy >= 1:
        #List array gives coefficients in the form b_0, b_1,...., b_NL,  Ra_0, Ia_0,..., Ra_NL, Ia_NL
        for r in range(8, 22):
            if r % 2 == 0:
                Real_coef.append(endf_dict[6][2]["subsection"][1]["A"][i][r])
            else:
                Imaginary_coef.append(endf_dict[6][2]["subsection"][1]["A"][i][r])
        for r in range(1, 8):
            b_coef.append(endf_dict[6][2]["subsection"][1]["A"][i][r])
        eta_val = eta_eval(1, 1, Energy, 1)
        dict_trig_pv = Trig_dict_production(len(Imaginary_coef)+1, eta_val, 1)
        dict_trig_nv = Trig_dict_production(len(Imaginary_coef)+1, eta_val, -1)
        leg_dict = Legendre_dict_production(2*len(Imaginary_coef)+1)
        Flag_back = True
        Back_cut_off_rate = 0
        for ang in angle_discretize_pos:
            #Ignores cross section data with scattering angle below both cut offs
            if CM_to_Lab_Frame(ang, Energy, 1, 1)<= max(cut_off,multiple_scattering_sd(Energy, 1, 1)*Gaussian_SD_cut_off):
                pass
            #Once backscatter cutoff is hit, backscattering is no longer possible and both backscatter and forward scatter 
            #density is given to the forward scatter proton
            elif ang >= math.cos(back_scatter_cut_off_cm) and Flag_back:
                Flag_back = False
                Back_cut_off_rate=Total_rate[-1]
                Total_rate.append(max(Back_cut_off_rate+
                    2*(Total_Rate_calc(
                        1, 1, Energy, 1, 1, ang, 1, Real_coef, Imaginary_coef, b_coef, dict_trig_pv, dict_trig_nv, leg_dict
                    )-Back_cut_off_rate),0)
                )
                yy.append(ang)
            elif not Flag_back:
                Total_rate.append(max(Back_cut_off_rate+
                    2*(Total_Rate_calc(
                        1, 1, Energy, 1, 1, ang, 1, Real_coef, Imaginary_coef, b_coef, dict_trig_pv, dict_trig_nv, leg_dict
                    )-Back_cut_off_rate),0)
                )
                yy.append(ang)
            else:
                Total_rate.append(max(
                    Total_Rate_calc(
                        1, 1, Energy, 1, 1, ang, 1, Real_coef, Imaginary_coef, b_coef, dict_trig_pv, dict_trig_nv, leg_dict
                    ),0)
                )
                yy.append(ang)
        if len(Total_rate) == 0 or Total_rate[-1] == 0:
            yy = [0.0]
            Total_rate = [1e-10]
        ErrorCheck(Total_rate, 1)
        Energy_discrete.append(Energy)
        Total_rate_True = []
        angle_discrete_True = []
        Flag_back_signed = True
        for i in range(len(Total_rate)):
            #Once backscatter cutoff is hit, backscattering is no longer possible and both backscatter and forward scatter 
            #density is given to the forward scatter proton
            if Flag_back_signed and yy[i] >= math.cos(back_scatter_cut_off_cm):
                angle_discrete_True.append(CM_to_Lab_Frame(yy[i], Energy, 1, 1))
                Total_rate_True.append(Total_rate[i])
                Flag_back_signed = False
            elif not Flag_back_signed:
                angle_discrete_True.append(CM_to_Lab_Frame(yy[i], Energy, 1, 1))
                Total_rate_True.append(Total_rate[i])
            else :
                angle_discrete_True.append(CM_to_Lab_Frame(yy[i], Energy, 1, 1))
                angle_discrete_True.insert(0,CM_to_Lab_Frame(-yy[i], Energy, 1, 1))
                Total_rate_True.insert(0, -Total_rate[i])
                Total_rate_True.append(Total_rate[i]) 
        if len(Total_rate_True) == 0 or Total_rate_True[-1] == 0:
            angle_discrete_True = [0.0]
            Total_rate_True = [1e-10]
        else:
            #normalizing constant
            max_rate = -Total_rate_True[0]
            for i in range(len(Total_rate_True)):
                Total_rate_True[i] += max_rate
        density[Energy] = [angle_discrete_True, Total_rate_True]
for key in density.keys():
    if len(density[key][0]) != len(density[key][1]):
     print(len(density[key][0]), len(density[key][1]))
f = open("Splines/hydrogen_el_ruth_cross_sec.txt", "w")
tmp = " ".join([str(z) for z in Energy_discrete])
f.write(tmp + "\n")
bool = True
for key in density.keys():
    tmp = " ".join([str(z) for z in density[key][0]])
    f.write(tmp + "\n")
    tmp = " ".join([str(z) for z in density[key][1]])
    f.write(tmp + "\n")
    bool = False
f.close()
