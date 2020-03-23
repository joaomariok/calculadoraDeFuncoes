#include <iostream>
#include <math.h>

#define THRESHOLD 1e-6
#define PRECISION 1e-9


using namespace std;


// Classe para funcoes do primeiro grau
class FuncaoDePrimeiroGrau {
public:
    double interseccaoX;
    double interseccaoY;
    bool ehCrescente;

    FuncaoDePrimeiroGrau(double a, double b)
    {
        this->a = a;
        this->b = b;
        ehCrescente = (a > 0) ? true : false;
        calcularInterseccoes();
    }

    void calcularInterseccoes()
    {
        interseccaoX = -b/a;
        interseccaoY = b;
    }

    void mostrarResultado()
    {
        cout << "f(x) = (" << a << ")x + (" << b << ")\n" << endl;

        if ( ehCrescente )
            cout << "FUNCAO CRESCENTE ( / )" << endl;
        else
            cout << "FUNCAO DECRESCENTE ( \\ )" << endl;

        cout << "INTERSECCAO COM X (y=0) : x = " << interseccaoX << endl;
        cout << "INTERSECCAO COM Y (x=0) : y = " << interseccaoY;
    }

private:
    double a;
    double b;
};


// Classe para funcoes do segundo grau
class FuncaoDeSegundoGrau {

public:
    double delta;
    double raizes[2];
    double extremoX;
    double extremoY;
    double derivada[2];
    bool raizesReais;
    bool raizesIguais;
    bool ehConcava;

    FuncaoDeSegundoGrau(double a, double b, double c)
    {
        this->a = a;
        this->b = b;
        this->c = c;
        calcularDelta();
        calcularRaizes();
        calcularPontoExtremo();
        calcularDerivada();
    }

    double calcularDelta()
    {
        delta = (pow(b,2) - (4*a*c));
        return delta;
    }

    void calcularRaizes()
    {
        if (delta < -THRESHOLD)
            raizesReais = false;
        else
        {
            if (fabs(delta) < THRESHOLD)
                delta = 0;

            raizesReais = true;
            raizes[0] = (-b + sqrt(delta))/(2*a);
            raizes[1] = (-b - sqrt(delta))/(2*a);
        }

        raizesIguais = (raizes[0] == raizes[1]) ? true : false;
    }

    void calcularPontoExtremo()
    {
        ehConcava = (a < 0) ? false : true;
        extremoX = -b/(2*a);
        extremoY = -delta/(4*a);
    }

    void calcularDerivada()
    {
        derivada[0] = b;
        derivada[1] = 2*a;
    }

    void mostrarResultado()
    {
        cout << "f(x)  = (" << a << ")x^2 + (" << b << ")x + (" << c << ")\n" << endl;
        cout << "f'(x) = (" << derivada[1] << ")x + (" << derivada[0] << ")\n" << endl;
        if ( !raizesReais ) // sem raizes reais
        {
            cout << "SEM RAIZES REAIS";
        }
        else // raizes reais
        {
            if ( raizesIguais ) // mesmas raizes
                cout << "RAIZ DUPLA:\n\n  x = " << raizes[0];
            else // raizes reais diferentes
                cout << "RAIZES:\n\n  x1 = " << raizes[0] << "\n  x2 = " << raizes[1];

            if( ehConcava )
                cout << "\n\n  Ponto de Minimo = (" << extremoX << ", " << extremoY << ")";
            else
                cout << "\n\n  Ponto de Maximo = (" << extremoX << ", " << extremoY << ")";
        }
    }

private:
    double a;
    double b;
    double c;
};


// Classe para funcoes do terceiro grau
class FuncaoDeTerceiroGrau {
public:
    double raizes[3];
    double derivada[3];
    double delta;
    bool raizesReais;
    bool raizesIguais;

    FuncaoDeTerceiroGrau(double a, double b, double c, double d)
    {
        this->a = a;
        this->b = b;
        this->c = c;
        this->d = d;
        calcularDelta();
        calcularDerivada();
        calcularRaizes();
    }

    void calcularDerivada()
    {
        derivada[0] = c;
        derivada[1] = 2*b;
        derivada[2] = 3*a;
    }

    void calcularPrimeiraRaiz()
    {
        double xAtual = 0;
        double xAnterior = -1;
        while (fabs(xAtual - xAnterior) > PRECISION) {
            xAnterior = xAtual;
            xAtual = xAnterior - (a*pow(xAnterior,3) + b*pow(xAnterior,2) + c*xAnterior + d)/
                                    (derivada[2]*pow(xAnterior,2) + derivada[1]*xAnterior + derivada[0]);
        }
        raizes[0] = xAtual;
    }

    double calcularDelta()
    {
        delta = 18*a*b*c*d - 4*pow(b,3)*d + pow(b,2)*pow(c,2) - 4*a*pow(c,3) - 27*pow(a,2)*pow(d,2);
        return delta;
    }

    void calcularRaizes()
    {
        calcularPrimeiraRaiz();

        if ( delta < -THRESHOLD )
        {
            raizesReais = false;
            raizesIguais = false;
        }
        else
        {
            briotRuffini();

            if ( fabs(delta) < THRESHOLD )
            {
                delta = 0;
                raizesReais = true;
                raizesIguais = true;
            }
            else
            {
                raizesReais = true;
                raizesIguais = false;
            }
        }
    }

    void briotRuffini ()
    {
        double resto = 0;
        double coefEqSegundoGrau[3];
        double coefEqTerceiroGrau[4] = {d, c, b, a};

        coefEqSegundoGrau[2] = coefEqTerceiroGrau[3];
        coefEqSegundoGrau[1] = raizes[0]*coefEqSegundoGrau[2] + coefEqTerceiroGrau[2];
        coefEqSegundoGrau[0] = raizes[0]*coefEqSegundoGrau[1] + coefEqTerceiroGrau[1];
        resto = raizes[0]*coefEqSegundoGrau[0] + coefEqTerceiroGrau[0];

        if (fabs(resto) > THRESHOLD)
            cout << "\n\nERRO!!! (sobrou resto) "<< resto <<"\n\n";

        FuncaoDeSegundoGrau funcaoAuxiliar(coefEqSegundoGrau[2], coefEqSegundoGrau[1], coefEqSegundoGrau[0]);
        raizes[1] = funcaoAuxiliar.raizes[0];
        raizes[2] = funcaoAuxiliar.raizes[1];
    }

    void mostrarResultado()
    {
        cout << "f(x)  = (" << a << ")x^3 + (" << b << ")x^2 + (" << c << ")x + (" << d << ")\n" << endl;
        cout << "f'(x) = (" << derivada[2] << ")x^2 + (" << derivada[1] << ")x + (" << derivada[0] << ")\n" << endl;
        cout << "RAIZES:" << endl << endl;
        if (raizesReais)
        {
            cout << "  x1 = " << raizes[0] << "\n  x2 = " << raizes[1] << "\n  x3 = " << raizes[2] << endl;
        }
        else
            cout << "x = " << raizes[0] << " e mais 2 raizes nao reais" << endl;
    }

private:
    double a;
    double b;
    double c;
    double d;
};


// Funcao main
int main()
{
    char opcao;

    cout << "<< CALCULADORA DE FUNCOES >>\n";

    while(1)
    {
        double a;
        double b;
        double c;
        double d;
        cout << "===============" << endl;
        cout << "Qual o grau da sua funcao (1 a 3)?" << endl;
        cout << "(Digite < 0 > para finalizar)" << endl;
        cout << "Opcao: ";
        cin >> opcao;

        if(opcao == '0')
            break;

        switch(opcao)
        {
        case '1':
            {
                cout << "\n\nENTRE COM OS VALORES < a > e < b >:\n\n";
                cout << "  f(x) = < a >*x + < b >\n\n";
                cout << "===============\n\na = ";
                cin >> a;

                if(a == 0)
                    break;

                cout << "b = ";
                cin >> b;
                cout << "\n---------------\n" << endl;

                FuncaoDePrimeiroGrau minhaFuncao1(a, b);
                minhaFuncao1.mostrarResultado();

                cout << "\n\n";
                break;
            }

        case '2':
            {
                cout << "\n\nENTRE COM OS VALORES < a >, < b > e < c >:\n\n";
                cout << "  f(x) = < a >*x^2 + < b >*x + < c >\n\n";
                cout << "===============\n\na = ";
                cin >> a;

                if(a == 0)
                    break;

                cout << "b = ";
                cin >> b;
                cout << "c = ";
                cin >> c;
                cout << "\n---------------\n" << endl;

                FuncaoDeSegundoGrau minhaFuncao2(a, b, c);
                minhaFuncao2.mostrarResultado();

                cout << "\n\n";
                break;
            }

        case '3':
            {
                cout << "\n\nENTRE COM OS VALORES < a >, < b >, < c > e < d >:\n\n";
                cout << "  f(x) = < a >*x^3 + < b >*x^2 + < c >*x + < d >\n\n";
                cout << "===============\n\na = ";
                cin >> a;

                if(a == 0)
                    break;

                cout << "b = ";
                cin >> b;
                cout << "c = ";
                cin >> c;
                cout << "d = ";
                cin >> d;
                cout << "\n---------------\n" << endl;

                FuncaoDeTerceiroGrau minhaFuncao3(a, b, c, d);
                minhaFuncao3.mostrarResultado();

                cout << endl;
                break;
            }

        default:
            {
                cout << "\n\nOPCAO INVALIDA!\n\n";
                break;
            }
        }
    }

    return 0;
}

