//
//  TMRSDarcyMemory.cpp
//  LinearTracer
//
//  Created by Omar Durán on 10/10/19.
//

#include "TMRSDarcyMemory.h"

TMRSDarcyMemory::TMRSDarcyMemory() {
    
    m_kappa.Resize(3,3);
    m_kappa.Identity();
    m_kappa_inv.Resize(3,3);
    m_kappa_inv.Identity();
    m_phi = 1.0;
    m_p = 0.0;
    m_p_n = 0.0;
    m_flux.Resize(3, 0.0);
    m_d = 0.0;
    
}

TMRSDarcyMemory::TMRSDarcyMemory(const TMRSDarcyMemory & other){
    m_kappa     =   other.m_kappa;
    m_kappa_inv = other.m_kappa_inv;
    m_phi       = other.m_phi;
    m_p         = other.m_p;
    m_p_n       = other.m_p_n;
    m_flux      = other.m_flux;
    m_d         = other.m_d;
}

const TMRSDarcyMemory & TMRSDarcyMemory::operator=(const TMRSDarcyMemory & other){
    // check for self-assignment
    if(&other == this){
        return *this;
    }
    m_kappa     =   other.m_kappa;
    m_kappa_inv = other.m_kappa_inv;
    m_phi       = other.m_phi;
    m_p         = other.m_p;
    m_p_n       = other.m_p_n;
    m_flux      = other.m_flux;
    m_d         = other.m_d;
    return *this;
}

TMRSDarcyMemory::~TMRSDarcyMemory(){

}

const std::string TMRSDarcyMemory::Name() const{
    return "TMRSDarcyMemory";
}

void TMRSDarcyMemory::Write(TPZStream &buf, int withclassid) const{
    buf.Write(m_kappa);
    buf.Write(m_kappa_inv);
    buf.Write(&m_phi);
    buf.Write(&m_p);
    buf.Write(&m_p_n);
    buf.Write(m_flux);
    buf.Write(&m_d);
}

void TMRSDarcyMemory::Read(TPZStream &buf, void *context){
    m_kappa.Read(buf,0);
    m_kappa_inv.Read(buf,0);
    buf.Read(&m_phi);
    buf.Read(&m_p);
    buf.Read(&m_p_n);
    buf.Read(m_flux);
    buf.Read(&m_d);
}

void TMRSDarcyMemory::Print(std::ostream &out) const{
    out << Name();
    m_kappa.Print(out);
    m_kappa_inv.Print(out);
    out << "\n m_phi = " << m_phi;
    out << "\n m_p = " << m_p;
    out << "\n m_p_n = " << m_p_n;
    out << "\n m_flux = " << m_flux;
    out << "\n m_d = " << m_d;
}

int TMRSDarcyMemory::ClassId() const{
    return Hash("TMRSDarcyMemory");
}
