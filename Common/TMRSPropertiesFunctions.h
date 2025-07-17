//
//  TMRSPropertiesFunctions.h
//  Buckley_Levertt
//
//  Created by Omar Durán on 6/16/20.
//

#ifndef TMRSPropertiesFunctions_h
#define TMRSPropertiesFunctions_h

#include "pzreal.h"
#include "TRMSpatialPropertiesMap.h"
#include <random>

class TMRSPropertiesFunctions
{
    public:
    
    /// Enumerate defining the function type
    enum EFunctionType { EConstantFunction = 0, EPiecewiseFunction = 1, ECircleLevelSetFunction = 2, EUNISIMFunction = 3, ESPECase10Function = 4, ERandomCirclesFunction = 5 };


    TMRSPropertiesFunctions(){
        m_function_type_kappa = EConstantFunction;
        m_function_type_phi = EConstantFunction;
        m_function_type_s0 = EConstantFunction;
        m_constant_val = 0.0;
    }
    TMRSPropertiesFunctions &operator=(TMRSPropertiesFunctions &other)
    {
        m_function_type_kappa = other.m_function_type_kappa;
        m_function_type_phi = other.m_function_type_phi;
        m_function_type_s0 = other.m_function_type_s0;
        return *this;
    };
    TMRSPropertiesFunctions(TMRSPropertiesFunctions &other)
    {
        m_function_type_kappa = other.m_function_type_kappa;
        m_function_type_phi = other.m_function_type_phi;
        m_function_type_s0 = other.m_function_type_s0;
    };
    ~TMRSPropertiesFunctions(){
        
    }
    
    void set_function_type_kappa(EFunctionType function_type_kappa){
        m_function_type_kappa   = function_type_kappa;
    }
    
    void set_function_type_phi(EFunctionType function_type_phi){
        m_function_type_phi     = function_type_phi;
    }
    
    /// @brief 
    /// @param function_type_s0 - function identifier
    /// @param value - optional to set a constant value if EConstantFunction is selected
    void set_function_type_s0(EFunctionType function_type_s0, REAL value = 0.0){
        m_function_type_s0      = function_type_s0;
        m_constant_val = value;
    }
    
    std::function<std::vector<REAL>(const TPZVec<REAL> & )> Create_Kappa_Phi(TRMSpatialPropertiesMap & map){
        
        return [& map] (const TPZVec<REAL> & pt) -> std::vector<REAL> {
            std::vector<REAL> kappa_and_phi;
            TPZManVector<REAL,3> x(pt);
            map.SampleKappaAndPhi(x,kappa_and_phi);
            return kappa_and_phi;
        };
        
    }
    
    std::function<std::vector<REAL>(const TPZVec<REAL> & )> Create_Kappa_Phi(){
        
        return [] (const TPZVec<REAL> & pt) -> std::vector<REAL> {
            std::vector<REAL> kappa_and_phi;
            REAL k_c,phi_c;
            k_c = 1.0e-7;
            phi_c = 0.1;
            REAL x = pt[0];
            REAL y = pt[1];
            REAL z = pt[2];
            REAL kx = k_c ;//* fabs(std::cos(0.2*x)*std::sin(0.1*y)*std::sin(0.1*z)) + k_c;
            REAL ky = 10.0 * k_c ;//* fabs(std::sin(0.1*x)*std::cos(0.2*y)*std::sin(0.1*z)) + k_c;
            REAL kz = 10.0 * k_c ;//* fabs(std::sin(0.1*x)*std::sin(0.1*y)*std::cos(0.2*z)) + k_c;
            REAL phi = phi_c ;//* fabs(std::cos(0.2*x)*std::cos(0.3*y)*std::cos(0.1*z)) + phi_c;
            kappa_and_phi.push_back(kx);
            kappa_and_phi.push_back(ky);
            kappa_and_phi.push_back(kz);
            kappa_and_phi.push_back(phi);
            return kappa_and_phi;
        };
        
    }

    std::function<REAL(const TPZVec<REAL> & )> Create_Kx(){
        
        switch (m_function_type_kappa) {
            case EConstantFunction:
                {
                    return [] (const TPZVec<REAL> & pt) -> REAL {
                            REAL x,y,kx;
                            x = pt[0];
                            y = pt[1];
                            kx = 1.0e-7;
                            return kx;
                        };
                }
                break;
            case ECircleLevelSetFunction:
                {
                    return [] (const TPZVec<REAL> & pt) -> REAL {
                            REAL x,y,kx,r,c,f;
                            r = 0.25;
                            c = 10;
                            x = pt[0];
                            y = pt[1];
                            f = -r*r + (-5 + x)*(-5 + x) + (-5 + y)*(-5 + y) - c;
                            if(f<0){
                                kx = 1.0e-10;
                            }else{
                                kx = 1.0e-7;
                            }
                            return kx;
                        };
                }
                break;
            case EPiecewiseFunction:
            {
                return [] (const TPZVec<REAL> & pt) -> REAL {
                        REAL y,kx;
                        y = pt[1];
                        if(y>5){
                             kx = 1.0e-10;
                         }else{
                             kx = 1.0e-7;
                         }
                         return kx;
                    };
            }
            break;
            default:
            {
                std::cout << " Function not implemented " << std::endl;
                DebugStop();
                return [] (const TPZVec<REAL> & pt) -> REAL {
                         return 0;
                     };
            }
                break;
        }
        
    }

    std::function<REAL(const TPZVec<REAL> & )> Create_Ky(){
        
        switch (m_function_type_kappa) {
            case EConstantFunction:
                {
                    return [] (const TPZVec<REAL> & pt) -> REAL {
                            REAL x,y,ky;
                            x = pt[0];
                            y = pt[0];
                            ky = 1.0e-7;
                            return ky ;
                        };
                }
                break;
            case ECircleLevelSetFunction:
                {
                    return [] (const TPZVec<REAL> & pt) -> REAL {
                            REAL x,y,ky,r,c,f;
                            r = 0.25;
                            c = 10;
                            x = pt[0];
                            y = pt[1];
                            f = -r*r + (-5 + x)*(-5 + x) + (-5 + y)*(-5 + y) - c;
                            if(f<0){
                                 ky = 1.0e-10;
                             }else{
                                 ky = 1.0e-7;
                             }
                             return ky;
                        };
                }
                break;
                case EPiecewiseFunction:
                {
                    return [] (const TPZVec<REAL> & pt) -> REAL {
                            REAL y,ky;
                            y = pt[1];
                            if(y>5){
                                 ky = 1.0e-10;
                             }else{
                                 ky = 1.0e-7;
                             }
                             return ky;
                        };
                }
                break;
            default:
            {
                std::cout << " Function not implemented " << std::endl;
                DebugStop();
                return [] (const TPZVec<REAL> & pt) -> REAL {
                         return 0;
                     };
            }
                break;
        }
        
    }
    
    std::function<REAL(const TPZVec<REAL> & )> Create_Kz(){
        
        switch (m_function_type_kappa) {
            case EConstantFunction:
                {
                    return [] (const TPZVec<REAL> & pt) -> REAL {
                            REAL x,y,kz;
                            x = pt[0];
                            y = pt[1];
                            kz = 1.0e-6;
                            return kz ;
                        };
                }
                break;
            case ECircleLevelSetFunction:
                {
                    return [] (const TPZVec<REAL> & pt) -> REAL {
                            REAL x,y,kz,r,c,f;
                            r = 0.25;
                            c = 10;
                            x = pt[0];
                            y = pt[1];
                            f = -r*r + (-5 + x)*(-5 + x) + (-5 + y)*(-5 + y) - c;
                            if(f<0){
                                 kz = 1.0e-10;
                             }else{
                                 kz = 1.0e-7;
                             }
                             return kz;
                        };
                }
                break;
                case EPiecewiseFunction:
                {
                    return [] (const TPZVec<REAL> & pt) -> REAL {
                            REAL y,kz;
                            y = pt[1];
                            if(y>5){
                                 kz = 1.0e-10;
                             }else{
                                 kz = 1.0e-7;
                             }
                             return kz;
                        };
                }
                break;
            default:
            {
                std::cout << " Function not implemented " << std::endl;
                DebugStop();
                return [] (const TPZVec<REAL> & pt) -> REAL {
                         return 0;
                     };
            }
                break;
        }
        
    }
    
    std::function<REAL(const TPZVec<REAL> & )> Create_phi(){
        
        switch (m_function_type_phi) {
            case EConstantFunction:
                {
                    return [] (const TPZVec<REAL> & pt) -> REAL {
                            REAL x,y,phi;
                            x = pt[0];
                            y = pt[1];
                            phi = 0.1;
                            return phi;
                        };
                }
                break;
            case ECircleLevelSetFunction:
                {
                    return [] (const TPZVec<REAL> & pt) -> REAL {
                            REAL x,y,phi,r,c,f;
                            r = 0.25;
                            c = 10;
                            x = pt[0];
                            y = pt[1];
                            f = -r*r + (-5 + x)*(-5 + x) + (-5 + y)*(-5 + y) - c;
                            if(f<0){
                                phi = 0.01;
                            }else{
                                phi = 0.1;
                            }
                            return phi;
                        };
                }
                break;
                case EPiecewiseFunction:
                {
                    return [] (const TPZVec<REAL> & pt) -> REAL {
                            REAL y,phi;
                            y = pt[1];
                            if(y<5){
                                 phi = 0.01;
                             }else{
                                 phi = 0.1;
                             }
                             return phi;
                        };
                }
                break;
            default:
            {
                std::cout << " Function not implemented " << std::endl;
                DebugStop();
                return [] (const TPZVec<REAL> & pt) -> REAL {
                         return 0;
                     };
            }
                break;
        }
        
    }
    
    std::function<REAL(const TPZVec<REAL> & )> Create_s0(){
        
        switch (m_function_type_s0) {
            case EConstantFunction:
                {
                    REAL s0 = m_constant_val;
                    return [s0] (const TPZVec<REAL> & pt) -> REAL {
                            REAL x,y,s;
                            x = pt[0];
                            y = pt[1];
                            s = s0;
                            return s;
                        };
                }
                break;
            case ECircleLevelSetFunction:
                {
                    return [] (const TPZVec<REAL> & pt) -> REAL {
                            REAL x,y,s,r,c,f;
                            r = 0.15; // mm
                            c = 0;
                            x = pt[0];
                            y = pt[1];
                            f = -r*r + (-202.65 + x)*(-202.65 + x) + (-800 + y)*(-800 + y) - c;
                            if(f<0){
                                s = 0.0;
                            }else{
                                s = 1.0;
                            }
                            return s;
                        };
                }
                break;
            case EPiecewiseFunction:
                {
                    return [] (const TPZVec<REAL> & pt) -> REAL {
                            REAL y,s;
                            y = pt[1];
                            if(y<0.3){
                                s = 1.0;
                            }else if (y>=0.3 && y<0.7){
                                s = 0.0;
                            }else{
                                s = 1.0;
                             }
                             return s;
                        };
                }
                break;
            case ERandomCirclesFunction:
                {
                    return [](const TPZVec<REAL> &pt) -> REAL
                    {
                        REAL x, y, s, c, f;
                        TPZVec<REAL> radius = {0.10, 0.05, 0.10, 0.05, 0.10}; // mm

                        c = 0;
                        TPZFNMatrix<100, REAL> center(5, 2);
                        center(0, 0) = 0.50;    center(0, 1) = 0.35;
                        center(1, 0) = 0.15;    center(1, 1) = 0.65;
                        center(2, 0) = 0.20;    center(2, 1) = 0.2;
                        center(3, 0) = 0.75;    center(3, 1) = 0.3;
                        center(4, 0) = 0.80;    center(4, 1) = 0.6;
                        x = pt[0];
                        y = pt[1];

                        for (int i = 0; i < center.Rows(); i++)
                        {
                            s = 1.0;
                            REAL r = radius[i];
                            f = -r * r + (x - center(i, 0)) * (x - center(i, 0)) + (y - center(i, 1)) * (y - center(i, 1)) - c;
                            if (f <= 0.0)
                            {
                                s = 0.0;
                                break; // If one circle is found, we can stop checking
                            }
                        }
                        return s;
                    };
                }
                break;
            default:
            {
                std::cout << " Function not implemented " << std::endl;
                DebugStop();
                return [] (const TPZVec<REAL> & pt) -> REAL {
                         return 0;
                     };
            }
                break;
        }
        
    }
    
    private:
    
    EFunctionType m_function_type_kappa;
    EFunctionType m_function_type_phi;
    EFunctionType m_function_type_s0;
    REAL m_constant_val;
};

#endif /* TMRSPropertiesFunctions_h */
