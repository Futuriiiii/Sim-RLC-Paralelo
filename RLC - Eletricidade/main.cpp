#include <cmath>
#include <iostream>

using namespace std;

#ifndef COMPONENTES_H
#define COMPONENTES_H

class componentes{

    public:
        componentes();
        ~componentes();
        void setR();
        void setL();
        void setC();
        void set_iL();
        void set_vC();
        void set_iC();
        void set_All();
        long double getR();
        long double getL();
        long double getC();
        long double get_iL();
        long double get_vC();
        long double get_iC();
        long double unidade(long double valor);
        void retorna_unidade(long double valor);

    private:
        int tipo;
        long double R, L, C, iL, iC, vC;
};

#endif

componentes::componentes(){

}

componentes::~componentes(){

}

void componentes::setR(){
    cout << "Insira o valor do Resistor: ";
    cin >> R;
    cout << "Insira a unidade desejada: " << endl;
    R = unidade(R);
}

void componentes::setL(){
    cout << "Insira o valor do Indutor: ";
    cin >> L;
    cout << "Insira a unidade desejada: " << endl;
    L = unidade(L);
}

void componentes::setC(){
    cout << "Insira o valor do Capacitor: ";
    cin >> C;
    cout << "Insira a unidade desejada: " << endl;
    C = unidade(C);
}

void componentes::set_iL(){
    cout << "Insira o valor da corrente inicial do indutor: ";
    cin >> iL;
    cout << "Insira a unidade desejada: " << endl;
    iL = unidade(iL);
}

void componentes::set_vC(){
    cout << "Insira o valor da tensao inicial do capacitor: ";
    cin >> vC;
    cout << "Insira a unidade desejada: " << endl;
    vC = unidade(vC);
}

void componentes::set_iC(){
    iC  = iL - (vC/R);
}

void componentes::set_All(){
    setR();
    setL();
    setC();
    set_iL();
    set_vC();
    set_iC();
}

long double componentes::getR(){
    return R;
}

long double componentes::getL(){
    return L;
}

long double componentes::getC(){
    return C;
}

long double componentes::get_iL(){
    return iL;
}

long double componentes::get_vC(){
    return vC;
}

long double componentes::get_iC(){
    return iC;
}

long double componentes::unidade(long double valor){
    int n;
    long double unidade;

    cout << "1 - 10 ^12 - Tera" << endl;
    cout << "2 - 10 ^9 - Giga" << endl;
    cout << "3 - 10 ^6 - Mega" << endl;
    cout << "4 - 10 ^3 - Kilo" << endl;
    cout << "5 - 10 ^0 - Sem submultiplo" << endl;
    cout << "6 - 10 ^-3 - Mili" << endl;
    cout << "7 - 10 ^-6 - Micro" << endl;
    cout << "8 - 10 ^-9 - Nano" << endl;
    cout << "9 - 10 ^-12 - Pico" << endl;

    cout << "Opcao desejada: " << endl;
    cin >> n;

    switch(n) {
    case 1:
        unidade = pow(10, 12);
        break;

    case 2:
        unidade = pow(10, 9);
        break;

    case 3:
        unidade = pow(10, 6);
        break;

    case 4:
        unidade = pow(10, 3);
        break;

    case 5:
        unidade = pow(10, 0);
        break;

    case 6:

        unidade = pow(10, -3);
        break;

    case 7:
        unidade = pow(10, -6);
        break;

    case 8:
        unidade = pow(10, -9);
        break;

    case 9:
        unidade = pow(10, -12);
        break;

    default:
        break;
    }

    return valor = valor * unidade;
}

void componentes::retorna_unidade(long double valor){
    int unidade = 0;
    long double aux;

    if(valor >= 0){
        aux = valor;
    }else{
        aux = -valor;
    }

    while(1){
        if(aux * pow(10,3) < 1000){
            unidade++;
            aux = aux * pow(10,3);
        }else{
            break;
        }

    }
    switch (unidade)
    {
    case 0:
        cout << valor * pow(10,(unidade*3));
        break;

    case 1:
        cout << valor * pow(10,(unidade*3)) << "(m";
        break;

    case 2:
        cout << valor * pow(10,(unidade*3)) << "(u";
        break;

    case 3:
        cout << valor * pow(10,(unidade*3)) << "(n";
        break;

    default:
        break;
    }

}

int main(){
    componentes comp;
    long double R, L, C, iL, iC, vC;
    long double tm, vtm;
    float sigma, omega0, omegaD;
    // M_E -> euler

    // fazendo a analise nodal do circuito, obtemos que iC = -iR - iL

    // vai ler as cargas iniciais do capacitor e indutor e os valores dos componentes
    comp.set_All();

    // vai atribuir os valores lidos às respectivivas variáveis
    R = comp.getR();
    L = comp.getL();
    C = comp.getC();
    iL = comp.get_iL();
    vC = comp.get_vC();
    iC = comp.get_iC();

    // calcula sigma e omega0
    sigma = 1/(2*R*C);
    omega0 = 1/sqrt(L*C);
    cout << "Sigma = " << sigma << "(s^-1)" << endl;
    cout << "Omega0 = " << omega0 << "(rad/s)" << endl;

    // compara os valores de sigma e omega para saber o tipo do circuito

    if(sigma > omega0){ // Circuito Superamortecido
        cout << "CIRCUITO SUPERAMORTECIDO" << endl;

        long double s1 = 0, s2 = 0, A1 = 0, A2 = 0;

        // calcula s1 e s2
        s1 = -sigma + (sqrt(pow(sigma, 2) - pow(omega0, 2)));
        s2 = -sigma - (sqrt(pow(sigma, 2) - pow(omega0, 2)));

/*
        // para calcular A1 e A2 temos que:
        // (-A1*s1) + (-A2*s2) = iC/C e A1 + A2 = vC, isolando A2 na segunda equação e substituindo na primeira, temos:
        // A2 = vC - A1 -> -A1*s1 - vC*s2 + A1*s2 = iC/C, isolando A1, temos:
        // -A1*s1 + A1*s2 = iC/C + vC*s2, colocando A1 em evidencia:
        // A1(-s1+s2) = iC/C + vC*s2, portanto,
        // A1 = (iC/C + vC*s2)/(-s1+s2)
        A1 = ((iC/C) + (vC*s2))/(-s1+s2);

        // usaremos A1 + A2 = vC para achar o valor de A2
        A2 = vC - A1;


        // para achar tm devemos derivar v(t) e igualar a 0
        // s1*A1*e^(s1*tm) + s2*A2*e^(s2*tm) = 0
        // s1*A1*e^(s1*tm) = -s2*A2*e^(s2*tm)
        // e^((s1-s2)*tm) = -s2*A2/s1*A1, passando ln dos dois lado temos:
        // (s1-s2)*tm = ln(-s2*A2) - ln(s1*A1), isolando o tm temos:
        // tm = (ln(-s2*A2) - ln(s1*A1))/(s1-s2)
        tm = (log(-s2*A2) - log(s1*A1))/(s1-s2);
*/

        long double eq = (-vC + (R *iL))/(R*C);
        A2 = ((-s1*vC) + eq)/(-s1+s2);
        A1 = vC - A2;

        tm = (log10(-(A2*s2)/(A1*s1))/log10(M_E)) / (s1-s2);

        vtm = (A1 * pow(M_E, s1*tm)) + (A2 * pow(M_E, s2*tm));

        cout << "s1 = " << s1 << endl;
        cout << "s2 = " << s2 << endl;
        cout << "A1 = " << A1 << endl;
        cout << "A2 = " << A2 << endl;
        //cout << "tm = " << tm << "(s)" << endl;
        cout << "tm = ";
        comp.retorna_unidade(tm);
        cout << "s)" << endl;
        //cout << "V(tm) = " << vtm << "(V)" << endl;
        cout << "V(tm) = ";
        comp.retorna_unidade(vtm);
        cout << "V)" << endl;

    }

    else if(sigma == omega0){ // Circuito Criticamente Amortecido
        cout << "CIRCUITO CRITICAMENTE AMORTECIDO" << endl;
        long double A1 = 0, A2 = 0, tm = 0;

        // para esse tipo de circuito temos que o A2 sera sempre igual ao vC
        A2 = vC;

        // para calcular o A1 iremos usar a equacao: iC/C = (dvc/dt),t=0 , assim obteremos
        // (iC/C) = A1 + (-sigma*A2)
        // elaborando essa equação com A1 e A2 chegamos que:
        // A1 = (iC/C) + (sigma*A2)
        A1 = (iC/C) + (sigma*A2);

        // para achar tm devemos derivar v(t) e igualar a 0
        // A1*e^(-sigma*tm) - sigma*e^(-sigma*tm) * (A1*tm + A2) = 0, passando o lado negativo para a outra parcela:
        // A1*e^(-sigma*tm) = sigma*e^(-sigma*tm) * (A1*tm + A2), simplificando, ficamos com:
        // A1 = sigma * (A1*tm + A2), então:
        // A1 = sigma*A1*tm + sigma*A2, passa a parcela com A2 para o outro lado e dps deixa tm isolado:
        // tm = (A1 - sigma*A2)/(sigma*A1)
        // tm = (A1 - sigma*A2)/(sigma*A1), podemos separar isso em:
        // tm = (A1/sigma*A1) - ((sigma*A2)/(sigma*A1)), simplificando, obtemos:
        // tm = (1-sigma) - (A2/A1)
        tm = (1/sigma) - (A2/A1);

        vtm = ((A1 * tm) + A2) * (pow(M_E,(-sigma * tm)));

        cout << "s = " << (-sigma) << endl;
        cout << "A1 = " << A1 << endl;
        cout << "A2 = " << A2 << endl;
        //cout << "tm = " << tm << "(s)" << endl;
        cout << "tm = ";
        comp.retorna_unidade(tm);
        cout << "s)" << endl;
        //cout << "V(tm) = " << vtm << "(V)" << endl;
        cout << "V(tm) = ";
        comp.retorna_unidade(vtm);
        cout << "V)" << endl;
    }

    else if(sigma < omega0){ // Circuito Subamortecido
        cout << "CIRCUITO SUBAMORTECIDO" << endl;
        float B1, B2;

        omegaD = sqrt(pow(omega0,2) - pow(sigma,2));

        // para esse tipo de circuito temos que o B1 sera sempre igual ao vC
        B1 = vC;

        // para calcular o B2 iremos usar a equacao: iC/C = (dvc/dt),t=0 , assim obteremos, de forma geral: iC/C = (sigma*B1) + (B2*omegaD)
        // elaborando essa equação com B1 e B2 chegamos que:
        // (iC/C)  (-sigma*B1) = (B2*omegaD) e então:
        // B2 = ((iC/C) - (sigma*B1))/omegaD
        B2 = ((iC/C) - (sigma*B1))/omegaD;

        // para achar tm devemos derivar v(t) e igualar a 0
        // -sigma*e^(-sigma*tm) (B1*cos(omegaD*tm) + B2*sen(omegaD*tm)) + e^(-sigma*tm) (-B1*omegaD*sen(omegaD*tm) + B2*omegaD*cos(omegaD*tm)) = 0, passando o lado negativo para a outra parcela:
        // e^(-sigma*tm) (-B1*omegaD*sen(omegaD*tm) + B2*omegaD*cos(omegaD*tm)) = sigma*e^(-sigma*tm) (B1*cos(omegaD*tm) + B2*sen(omegaD*tm)), simplificando, temos:
        // (-B1*omegaD*sen(omegaD*tm) + B2*omegaD*cos(omegaD*tm)) = sigma * (B1*cos(omegaD*tm) + B2*sen(omegaD*tm)), multiplicando o sigma no lado direito:
        // -B1*omegaD*sen(omegaD*tm) + B2*omegaD*cos(omegaD*tm) = B1*sigma*cos(omegaD*tm) + B2*sigma*sen(omegaD*tm), deixa cos e sen em lados separados:
        // B2*omegaD*cos(omegaD*tm) - B1*sigma*cos(omegaD*tm) = B2*sigma*sen(omegaD*tm) + B1*omegaD*sen(omegaD*tm), deixa os cos e sen em evidencia:
        // cos(omegaD*tm) * (B2*omegaD - B1*sigma) = sen(omegaD*tm) * (B2*sigma + B1*omegaD), para chegar na tg:
        // tg(omegaD*tm) = (B2*omegaD - B1*sigma) / (B2*sigma + B1*omegaD), passando o arctg temos:
        // omegaD*tm = arctg((B2*omegaD - B1*sigma) / (B2*sigma + B1*omegaD)), portanto:
        // tm = arctg((B2*omegaD - B1*sigma) / (B2*sigma + B1*omegaD))/omegaD
        tm = atan((B2*omegaD - B1*sigma) / (B2*sigma + B1*omegaD))/omegaD;
        if(tm < 0){
            tm = tm + ((M_PI/omegaD)/2);
        }

        vtm = (pow(M_E, (-sigma * tm))) * (B1 * cos(omegaD * tm) + (B2 * sin(omegaD * tm)));

        cout << "OmegaD = " << omegaD << " (rad/s)" << endl;
        cout << "B1 = " << B1 << endl;
        cout << "B2 = " << B2 << endl;
        //cout << "tm = " << tm << "(s)" << endl;
        cout << "tm = ";
        comp.retorna_unidade(tm);
        cout << "s)" << endl;
        //cout << "V(tm) = " << vtm << "(V)" << endl;
        cout << "V(tm) = ";
        comp.retorna_unidade(vtm);
        cout << "V)" << endl;
    }

    cout << "\n*******AUTORES*******\nCarlos Rafael Torres Miranda Novack - 20210066961\nArthur Henrique da Silva - 20200077587" << endl;

    return 0;
}
