#include <iostream> // Fluxo de I/O por monitor, teclado...
#include <fstream>  /* i/ofstream | fstream é um manipulador de fluxos de dados de arquivos de computador
especializado para o tipo de dado nativo char. Ele permite ler e escrever em modo de texto (utiliza-se os operadores de
deslocamento de bits, << e >>) ou binário (utiliza-se os métodos read e write para buffers de dado). Fluxo de Arquivos*/
//#include <cmath>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>  //getchar
#include <string>
#include <stdlib.h>     //para usar o rand()
#include "MersenneTwister.h"

using namespace std;

//Início do código que tem o objetivo de gerar cadeias de polímeros com persistência dada no arquivo parametros.txt.
int main()
{
    //Espaço para declaração das variáveis usadas neste código:
    float densidade,beads,persistencia,volume,r,temposicao,totalposicoes;
    int i_seed,n,a,b,c,i_cadeias,ncadeias,x,y,z,sv,armadilha,svr,ocupadas,aresta,tentativas;
    int xma,xme,yma,yme,zma,zme,xt,yt,zt,perauxme,perauxyma,perauxyme,perauxzma,perauxzme,perauxma,f;
    unsigned int caixa[80][80][80],rede[80][80][80], ctemp[10000]; //retira o espaço de armazenamento para valores negativos, aumentando o
    //espaço das matrizes caixa e rede.
	//ctemp tem relação direta com caixa, rede e densidade. Antes de qualquer alteração analise o que pode acontecer.
	
    MTRand mtrand1; // Algum tipo de declaração da MersenneTwister.h.

/*---------------------------------------------------------------------------------------------------------------------------
/ Apontando os arquivos de saída */

//Setando teste como objeto para o arquivo permanente.dat
ofstream permanente; // cria o objeto permanente.
permanente.open("caixamc.dat"); //abre o arquivo permanente.dat que deve estar no diretório corrente.

/*---------------------------------------------------------------------------------------------------------------------------*/
//Recuperação das propriedades do polímero: */

    ifstream parametros; // cria o objeto leitura.
    parametros.open("parametros.txt"); //abre o arquivo parametros.txt, que deve estar no diretório corrente.

    // Leitura das propriedades do polímero.
    
    parametros.width(3);
    parametros>>densidade;
    parametros.width(3);
    parametros>>beads;
    parametros.width(3);
    parametros>>persistencia;
    parametros.width(3);
    parametros>>i_seed;
    parametros.width(3);
    parametros>>aresta;
    
    parametros.close(); //fecha o arquivo setado pelo objeto desenvolvimento.

/*-----------------------------------------------------------------------------------------------------------------
 Área de inicialização de variáveis e afins: */

volume=(aresta*aresta*aresta);
totalposicoes=0.0;
mtrand1.seed(i_seed); //aponta a semente para a gereção dos números aleatórios usados no decorrer do código.

/*-----------------------------------------------------------------------------------------------------------------*/

i_cadeias=((densidade*volume)/beads); /* Numero de cadeias necessárias para preencher a densidade indicada no arquivo
                                       parametros.txt. A função floor foi usada para arredondar para baixo o valor
                                       do seu argumento.*/
                                       
//-----------------------------------------------------------------------------------------------------------------
//Visualização das propriedades e afins:

cout<<"Densidade: "<<densidade<<endl;
cout<<"Beads: "<<beads<<endl;
cout<<"Persistencia: "<<persistencia<<endl;
cout<<"Aresta: "<<aresta<<endl;
cout<<"Semente: "<<i_seed<<endl;
cout<<"Numero de cadeias necessárias: "<<i_cadeias<<endl;

//------------------------------------------------------------------------------------------------------------------
// Zerando os sitios da cadeia principal:

for(a=0.0;a<=(aresta-1.0);a++)
{
  for(b=0.0;b<=(aresta-1.0);b++)
  {
    for(c=0.0;c<=(aresta-1.0);c++)
    {
      caixa[a][b][c]=0.0;
    }
  }
}

//-----------------------------------------------------------------------------------------------------------------
// Zerando contador do número total de armadilhas.
armadilha=0.0;

//-----------------------------------------------------------------------------------------------------------------
// Geração de cadeias considerando persistência indicada no parametros.txt.
//-----------------------------------------------------------------------------------------------------------------

for(ncadeias=1;ncadeias<=i_cadeias;ncadeias++)
{    
    reset:
	temposicao=0.0;
	//Setando temp como objeto para o arquivo temp.dat
	ofstream temp; // cria o objeto temp.
	temp.open("temp.dat"); //abre o arquivo temp.dat que deve estar no diretório corrente.
	/* É necessário setar o objeto temp neste ponto pois o arquivo temp.dat terá o ojetivo de armazenar uma cadeia temporaria 
	 * do conjunto de cadeais considerado. Quando acionamos reset ou adição de um novo bead o temp.dat é fechado e "ofstream temp"
	 * limpa todo o arquivo temp.dat e posiciona o cursor para adição de mais dados neste arquivo na primeira coluna na primeira linha. */
    
//------------------------------------------------------------------------------------------------------------------
// Zerando a cadeia de monômeros auxiliar.

    for(a=0;a<=(aresta-1.0);a++)
    {
        for(b=0;b<=(aresta-1.0);b++)
        {
            for(c=0;c<=(aresta-1.0);c++)
            {
                rede[a][b][c]=0.0; // zerando todos os s�tios da rede
            }
        }
    }

//--------------------------------------------------------------------------------------------------------------------------------
// Sorteio da origem na caixa de MC

    do
    {
        //sorteio de numeros para origem entre li e ls.
        x = mtrand1.randDblExc( aresta ); 
        y = mtrand1.randDblExc( aresta );
        z = mtrand1.randDblExc( aresta );
    } while (caixa[x][y][z]!=0.0);
	
    temp.width(3);
    temp<<x;
    temp.width(3);
    temp<<y;
    temp.width(3);
    temp<<z;
    
    xt=x;
    yt=y;
    zt=z;
    
    rede[x][y][z]=1.0; // primeiro sítio ocupado.
    
    temposicao=temposicao+1.0;
    
//---------------------------------------------------------------------------------------------------------------------------------
// Sorteio da primeira persistência

    r =  mtrand1.rand();

    if (r<1.0/6.0)
    {
        perauxma=1.0;
        xme=yma=yme=zma=zme=0.0;
        perauxme=perauxyma=perauxyme=perauxzma=perauxzme=0.0;
    }

    if ((r>1.0/6.0) and (r<=2.0/6.0))
    {
        perauxme=1.0;
        xma=yma=yme=zma=zme=0.0;
        perauxma=perauxyma=perauxyme=perauxzma=perauxzme=0.0;
    }

    if ((r>2.0/6.0) and (r<=3.0/6.0))
    {
        perauxyma=1.0;
        xma=xme=yme=zma=zme=0.0;
        perauxma=perauxme=perauxyme=perauxzma=perauxzme=0.0;
    }

    if ((r>3.0/6.0) and (r<=4.0/6.0))
    {
        perauxyme=1.0;
        xma=xme=yma=zma=zme=0.0;
        perauxma=perauxme=perauxyma=perauxzma=perauxzme=0.0;
    }

    if ((r>4.0/6.0) and (r<=5.0/6.0))
    {
        perauxzma=1.0;
        xma=xme=yma=yme=zme=0.0;
        perauxma=perauxme=perauxyma=perauxyme=perauxzme=0.0;
    }

    if ((r>5.0/6.0) and (r<=1.0))
    {
        perauxzme=1.0;
        xma=xme=yma=yme=zma=0.0;
        perauxma=perauxme=perauxyma=perauxyme=perauxzma=0.0;
    }

//--------------------------------------------------------------------------------------------------------------------------------
// Geração das posições dos beads da cadeia de numero dado pelo "ncadeias".

    for(n=1;n<beads;n++)
    {
        
        tentativas=0;
        
        start:
        
        tentativas=tentativas+1;
        
        if (tentativas==200) goto reset;
        
        r =  mtrand1.rand();
         
        ocupadas=0.0;

//-------------------------------------------------------------------------------------------------------------------------------------------------------
// Soma e teste de condição de vizinhos para checar se a geração se armadilhou. Checagem das posições ocupadas para conferência.

        sv = caixa[z+1][y][x] + caixa[z-1][y][x] + caixa[z][y+1][x] + caixa[z][y-1][x] + caixa[z][y][x+1] + caixa[z][y][x-1]; // soma dos vizinhos na caixa principal.
        svr = rede[z+1][y][x] + rede[z-1][y][x] + rede[z][y+1][x] + rede[z][y-1][x] + rede[z][y][x+1] + rede[z][y][x-1]; //soma dos vizinhos na rede.
        		
	    if (caixa[z+1][y][x]==1.0)
	    {
            ocupadas=ocupadas+1.0;
	    }

	    if (caixa[z-1][y][x]==1.0)
	    {
            ocupadas=ocupadas+1.0;
	    }

	    if (caixa[z][y+1][x]==1.0)
	    {
            ocupadas=ocupadas+1.0;
	    }

	    if (caixa[z][y-1][x]==1.0)
	    {
            ocupadas=ocupadas+1.0;
	    }

	    if (caixa[z][y][x+1]==1.0)
	    {
	    	ocupadas=ocupadas+1.0;
	    }

	    if (caixa[z][y][x-1]==1.0)
	    {
            ocupadas=ocupadas+1.0;
	    }


	    if (rede[z+1][y][x]==1.0)
	    {
            ocupadas=ocupadas+1.0;
	    }

	    if (rede[z-1][y][x]==1.0)
	    {
            ocupadas=ocupadas+1.0;
	    }

	    if (rede[z][y+1][x]==1.0)
	    {
            ocupadas=ocupadas+1.0;
	    }

	    if (rede[z][y-1][x]==1.0)
	    {
	        ocupadas=ocupadas+1.0;
	    }

	    if (rede[z][y][x+1]==1.0)
	    {
            ocupadas=ocupadas+1.0;
	    }

	    if (rede[z][y][x-1]==1.0)
	    {
	        ocupadas=ocupadas+1.0;
	    }

	    if (ocupadas==6.0) //soma dos vizinhos da rede e da caixa.
	    {
            armadilha++;
            if (armadilha%30==1) 
			{
				cout<<"Armadilha: "<<armadilha<<". Cadeia: "<<ncadeias<<". Vizinhos: "<<sv+svr<<endl;       
			}
            temp.close();
            goto reset;
	    }

//-------------------------------------------------------------------------------------------------------------------------------------------------------

        if (perauxma>=1.0)
        {
            if (r<=persistencia)
            {
                xt = x + 1.0;     
                y=yt;
                z=zt;
                xma=1.0;
                xme=yma=yme=zma=zme=0.0;
                perauxme=perauxyma=perauxyme=perauxzma=perauxzme=0.0;
            }

            if (r>persistencia && r<= (persistencia*4.0+1.0)/5.0)
            {
                yt = y + 1.0;     
                x=xt;
                z=zt;
                yma=1.0;
                xma=xme=yme=zma=zme=0.0;
                perauxma=perauxme=perauxyme=perauxzma=perauxzme=0.0;
            }

            if (r > (persistencia*4.0+1.0)/5. && r <= (persistencia*3.0+2.0)/5.0 )
            {
                xt = x - 1.0;    
                y=yt;
                z=zt;
                xme=1.0;
                xma=yma=yme=zma=zme=0.0;
                perauxma=perauxyma=perauxyme=perauxzma=perauxzme=0.0;
            }

            if ( r > (persistencia*3.0+2.0)/5.0 && r<= (persistencia*2.0+3.0)/5.0)
            {
                yt = y - 1.0;       
                x=xt;
                z=zt;
                yme=1.0;
                xma=xme=yma=zma=zme=0.0;
                perauxme=perauxma=perauxyma=perauxzma=perauxzme=0.0;
            }

            if (r > (persistencia*2.0+3.0)/5.0 && r<= (persistencia+4.0)/5.0)
            {
                zt = z +1.0;     
                y=yt;
                x=xt;
                zma=1.0;
                xma=xme=yma=yme=zme=0.0;
                perauxme=perauxma=perauxyma=perauxyme=perauxzme=0.0;
            }

            if ( r > (persistencia+4.0)/5.0 && r<=1.0)
            {
                zt = z -1.0;     
                y=yt;
                x=xt;
                zme=1.0;
                xma=xme=yma=yme=zma=0.0;
                perauxme=perauxma=perauxyma=perauxyme=perauxzma=0.0;
            }
        }

        if (perauxme>=1.0)
        {
            if (r<=persistencia)
            {
                xt = x - 1.0;    
                y=yt;
                z=zt;
                xme=1.0;
                xma=yma=yme=zma=zme=0.0;
                perauxma=perauxyma=perauxyme=perauxzma=perauxzme=0.0;
            }

            if (r>persistencia && r<=(persistencia*4.0+1.0)/5.0)
            {
                yt = y +1.0;     
                x=xt;
                z=zt;
                yma=1.0;
                xma=xme=yme=zma=zme=0.0;
                perauxma=perauxme=perauxyme=perauxzma=perauxzme=0.0;
            }

            if (r > (persistencia*4.0+1.0)/5.0 && r <= (persistencia*3.0+2.0)/5.0 )
            {
                xt = x +1.0;     
                y=yt;
                z=zt;
                xma=1.0;
                xme=yma=yme=zma=zme=0.0;
                perauxme=perauxyma=perauxyme=perauxzma=perauxzme=0.0;
            }

            if ( r > (persistencia*3.0+2.0)/5.0 && r<= (persistencia*2.0+3.0)/5.0)
            {
                yt = y - 1.0;       
                x=xt;
                z=zt;
                yme=1.0;
                xma=xme=yma=zma=zme=0.0;
                perauxme=perauxma=perauxyma=perauxzma=perauxzme=0.0;
            }

            if (r > (persistencia*2.0+3.0)/5.0 && r<= (persistencia+4.0)/5.0)
            {
                zt = z +1.0;     
                y=yt;
                x=xt;
                zma=1.0;
                xma=xme=yma=yme=zme=0.0;
                perauxme=perauxma=perauxyma=perauxyme=perauxzme=0.0;
            }

            if ( r > (persistencia+4.0)/5.0 && r<=1.0)
             {
                zt = z -1.0;     
                y=yt;
                x=xt;
                zme=1.0;
                xma=xme=yma=yme=zma=0.0;
                perauxme=perauxma=perauxyma=perauxyme=perauxzma=0.0;
            }
        }

        if (perauxyma>=1.0)
        {
            if (r<=persistencia)
            {
                yt = y +1.0;     
                x=xt;
                z=zt;
                yma=1.0;
                xma=xme=yme=zma=zme=0.0;
                perauxma=perauxme=perauxyme=perauxzma=perauxzme=0;
            }
            if (r>persistencia && r<=(persistencia*4.0+1.0)/5.0)
            {
                xt = x +1.0;     
                y=yt;
                z=zt;
                xma=1.0;
                xme=yma=yme=zma=zme=0.0;
                perauxme=perauxyma=perauxyme=perauxzma=perauxzme=0.0;
            }

            if (r > (persistencia*4.0+1.0)/5.0 && r <= (persistencia*3.0+2.0)/5.0 )
            {
                xt = x - 1.0;    
                y=yt;
                z=zt;
                xme=1.0;
                xma=yma=yme=zma=zme=0.0;
                perauxma=perauxyma=perauxyme=perauxzma=perauxzme=0.0;
            }

            if ( r > (persistencia*3.0+2.0)/5.0 && r<= (persistencia*2.0+3.0)/5.0)
            {
                yt = y - 1.0;       
                x=xt;
                z=zt;
                yme=1.0;
                xma=xme=yma=zma=zme=0.0;
                perauxme=perauxma=perauxyma=perauxzma=perauxzme=0.0;
            }

            if (r > (persistencia*2.0+3.0)/5.0 && r<= (persistencia+4.0)/5.0)
            {
                zt = z +1.0;     
                y=yt;
                x=xt;
                zma=1.0;
                xma=xme=yma=yme=zme=0.0;
                perauxme=perauxma=perauxyma=perauxyme=perauxzme=0.0;
            }

            if ( r > (persistencia+4.0)/5.0 && r<=1.0)
            {
                zt = z -1.0;     
                y=yt;
                x=xt;
                zme=1.0;
                xma=xme=yma=yme=zma=0.0;
                perauxme=perauxma=perauxyma=perauxyme=perauxzma=0.0;
            }
        }

        if (perauxyme>=1.0)
        {
            if (r<=persistencia)
            {
                yt = y - 1.0;       
                x=xt;
                z=zt;
                yme=1.0;
                xma=xme=yma=zma=zme=0.0;
                perauxme=perauxma=perauxyma=perauxzma=perauxzme=0.0;
            }

            if (r>persistencia && r<=(persistencia*4.0+1.0)/5.0)
            {
                yt = y +1.0;     
                x=xt;
                z=zt;
                yma=1.0;
                xma=xme=yme=zma=zme=0.0;
                perauxma=perauxme=perauxyme=perauxzma=perauxzme=0.0;
            }

            if (r > (persistencia*4.0+1.0)/5.0 && r <= (persistencia*3.0+2.0)/5.0 )
            {
                xt = x - 1.0;   
                y=yt;
                z=zt;
                xme=1.0;
                xma=yma=yme=zma=zme=0.0;
                perauxma=perauxyma=perauxyme=perauxzma=perauxzme=0.0;
            }

            if ( r > (persistencia*3.0+2.0)/5.0 && r<= (persistencia*2.0+3.0)/5.0)
            {
                xt = x +1.0;     
                y=yt;
                z=zt;
                xma=1.0;
                xme=yma=yme=zma=zme=0.0;
                perauxme=perauxyma=perauxyme=perauxzma=perauxzme=0.0;
            }

            if (r > (persistencia*2.0+3.0)/5.0 && r<= (persistencia+4.0)/5.0)
            {
                zt = z +1.0;     
                y=yt;
                x=xt;
                zma=1.0;
                xma=xme=yma=yme=zme=0.0;
                perauxme=perauxma=perauxyma=perauxyme=perauxzme=0.0;
            }

            if ( r > (persistencia+4.0)/5.0 && r<=1.0)
            {
                zt = z -1.0;     
                y=yt;
                x=xt;
                zme=1.0;
                xma=xme=yma=yme=zma=0.0;
                perauxme=perauxma=perauxyma=perauxyme=perauxzma=0;
            }
        }

        if (perauxzma>=1.0)
        {
            if (r<=persistencia)
            {
                zt = z +1.0;     
                y=yt;
                x=xt;
                zma=1.0;
                xma=xme=yma=yme=zme=0.0;
                perauxme=perauxma=perauxyma=perauxyme=perauxzme=0.0;
            }

            if (r>persistencia && r<=(persistencia*4.0+1.0)/5.0)
            {
                yt = y +1.0;     
                x=xt;
                z=zt;
                yma=1.0;
                xma=xme=yme=zma=zme=0.0;
                perauxma=perauxme=perauxyme=perauxzma=perauxzme=0.0;
            }

            if (r > (persistencia*4.0+1.0)/5.0 && r <= (persistencia*3.0+2.0)/5.0 )
            {
                xt = x - 1.0;    
                y=yt;
                z=zt;
                xme=1.0;
                xma=yma=yme=zma=zme=0.0;
                perauxma=perauxyma=perauxyme=perauxzma=perauxzme=0.0;
            }

            if ( r > (persistencia*3.0+2.0)/5.0 && r<= (persistencia*2.0+3.0)/5.0)
            {
                yt = y - 1.0;       
                x=xt;
                z=zt;
                yme=1.0;
                xma=xme=yma=zma=zme=0.0;
                perauxme=perauxma=perauxyma=perauxzma=perauxzme=0.0;
            }

            if (r > (persistencia*2.0+3.0)/5.0 && r<= (persistencia+4.0)/5.0)
            {
                xt = x +1.0;     
                y=yt;
                z=zt;
                xma=1.0;
                xme=yma=yme=zma=zme=0.0;
                perauxme=perauxyma=perauxyme=perauxzma=perauxzme=0.0;
            }

            if ( r > (persistencia+4.0)/5.0 && r<=1.0)
            {
                zt = z -1.0;     
                y=yt;
                x=xt;
                zme=1.0;
                xma=xme=yma=yme=zma=0.0;
                perauxme=perauxma=perauxyma=perauxyme=perauxzma=0.0;
            }
        }

        if (perauxzme>=1.0)
        {
            if (r<=persistencia)
            {
                zt = z -1.0;     
                y=yt;
                x=xt;
                zme=1.0;
                xma=xme=yma=yme=zma=0.0;
                perauxme=perauxma=perauxyma=perauxyme=perauxzma=0.0;
            }

            if (r>persistencia && r<=(persistencia*4.0+1.0)/5.0)
            {
                yt = y +1.0;     
                x=xt;
                z=zt;
                yma=1.0;
                xma=xme=yme=zma=zme=0.0;
                perauxma=perauxme=perauxyme=perauxzma=perauxzme=0.0;
            }

            if (r > (persistencia*4.0+1.0)/5.0 && r <= (persistencia*3.0+2.0)/5.0 )
            {
                xt = x - 1.0;    
                y=yt;
                z=zt;
                xme=1.0;
                xma=yma=yme=zma=zme=0.0;
                perauxma=perauxyma=perauxyme=perauxzma=perauxzme=0.0;
            }

            if ( r > (persistencia*3.0+2.0)/5.0 && r<= (persistencia*2.0+3.0)/5.0)
            {
                yt = y - 1.0;       
                x=xt;
                z=zt;
                yme=1.0;
                xma=xme=yma=zma=zme=0.0;
                perauxme=perauxma=perauxyma=perauxzma=perauxzme=0.0;
            }

            if (r > (persistencia*2.0+3.0)/5.0 && r<= (persistencia+4.0)/5.0)
            {
                zt = z +1.0;     
                y=yt;
                x=xt;
                zma=1.0;
                xma=xme=yma=yme=zme=0.0;
                perauxme=perauxma=perauxyma=perauxyme=perauxzme=0.0;
            }

            if ( r > (persistencia+4.0)/5.0 && r<=1.0)
            {
                xt = x +1.0;     
                y=yt;
                z=zt;
                xma=1.0;
                xme=yma=yme=zma=zme=0.0;
                perauxme=perauxyma=perauxyme=perauxzma=perauxzme=0.0;
            }
        }
       
//-------------------------------------------------------------------------------------------------------------------------------------------------------
// Condições periódicas de contorno.

        if (xt<0.0) xt=aresta-1.0;      
        if (yt<0.0) yt=aresta-1.0;
        if (zt<0.0) zt=aresta-1.0;
        if (xt>(aresta-1.0)) xt=0.0;
        if (yt>(aresta-1.0)) yt=0.0;
        if (zt>(aresta-1.0)) zt=0.0;
               
//Fim das condições periódicas de contorno.
//------------------------------------------------------------------------------------------------------------------------------------------------------

// Verificação da rede ou da caixa principal ocupada. Se a posição futura estiver ocupada a posição anterior será selecionada.

		if( rede[xt][yt][zt]==1.0 or caixa[xt][yt][zt]==1.0) // rede ocupada.
        {

            if (xma==1.0)
            {
                perauxma=perauxma+1.0;
				xt=x;
            }

            if (xme==1.0)
            {
               perauxme=perauxme+1.0;
               xt=x;
            }

            if (yma==1.0)
            {
                perauxyma=perauxyma+1.0;
                yt=y;
            }

            if (yme==1.0)
            {
                perauxyme=perauxyme+1.0;
                yt=y;
            }

            if (zma==1.0)
            {
                perauxzma=perauxzma+1.0;
                zt=z;
            }

            if (zme==1.0)
            {
            	perauxzme=perauxzme+1.0;
                zt=z;
            }

            goto start;

        }
 
 //Fim do teste que analisa se a posição futura está ocupada.
//--------------------------------------------------------------------------------------------------------------------------------------------------------

//Adição do bead ao polímero
        
        if( (rede[xt][yt][zt]) == 0.0 and (caixa[xt][yt][zt]== 0.0))
        {
            			
            x=xt;
            y=yt;
            z=zt;
           
            rede[x][y][z]=1.0;
            
            if (xma==1.0) perauxma=perauxma+1.0;
            if (xme==1.0) perauxme=perauxme+1.0;
            if (yma==1.0) perauxyma=perauxyma+1.0;
            if (yme==1.0) perauxyme=perauxyme+1.0;
            if (zma==1.0) perauxzma=perauxzma+1.0;
            if (zme==1.0) perauxzme=perauxzme+1.0;
            
            temp.width(3);
            temp<<x;
            temp.width(3);
            temp<<y;
            temp.width(3);
            temp<<z;   
            temposicao=temposicao+1.0;
		}
		
//Fim da adição dos beads.
//------------------------------------------------------------------------------------------------------------------------------------------------------
		
//-------------------------------------------------------------------------------------------------------------------------------------------------------
    } // Fim da geração de cadeias do loop da variável n.
	totalposicoes=totalposicoes+temposicao;
	temp.close(); //fecha o arquivo temp.dat que recebeu os beads da cadeia considerada.
	ifstream gravatemp; // cria o objeto gravatemp.
    gravatemp.open("temp.dat"); //gravatemp abre o arquivo temp.dat, que acabou de ser fechado com uma cadeia a ser adicionada a caixa principal.
	f=0.0;
	
// grava a cadeia considerada que está no arquivo temp.dat na matriz ctemp.	
	while (gravatemp.eof()==0.0)
	{		
		gravatemp.width(3);
		gravatemp>>ctemp[f];
		f++;	
	}
	gravatemp.close();

// grava a matriz ctemp[a] no arquivo setado pelo objeto permanente.	
	for (a=0.0;a<f;a++)
	{	
		permanente.width(3);
		permanente<<ctemp[a];
	}
	permanente<<endl;
	
//-------------------------------------------------------------------------------------------------------------------------------------------------------

//Adicionando a cadeia ncadeias de polimeros a caixa de MC:

    for(a=0;a<=(aresta-1);a++)
    {
        for(b=0;b<=(aresta-1);b++)
        {
            for(c=0;c<=(aresta-1);c++)
            {

                if (rede[a][b][c]==1.0)
                {
                    caixa[a][b][c]=rede[a][b][c]; // adicionado cadeia de polimeros a caixa de MC.
                }

            }
        }
    }
//--------------------------------------------------------------------------------------------------------------------------------------------
	
	if (ncadeias%30==1) 
    {
    	cout<<"Cadeia: "<<ncadeias<<endl;       
    }
    	
}//-------------------------------------------------------------------------------------------------------------------------------------------
// Fim do for que controla o numero de cadeias.
//--------------------------------------------------------------------------------------------------------------------------------------------
 permanente.close();
 cout<<"Última cadeia: "<<ncadeias-1<<" . Armadilhas totais: "<<armadilha<<endl;
 cout<<"Porcentagem de posições ocupadas: "<<(totalposicoes/(aresta*aresta*aresta))*100.0<<endl;
 cout<<endl;
 cout<<"----"<<endl;
 cout<<" OK. "<<endl;
 cout<<"----";
 return 0;
}
