/*    ******************************************************************
*     ***    Este programa simula uma rede de polímeros.             ***
*     ***              									             *** 
* 	  ******                                                      ******
*     ***  Referencia base: Arquivo-fonte sim.cpp fornecido          ***
*     ***  por Ronaldo Junio e fornecido  por  Vitor Barbante        ***
*     ***  Pereira Leite.                                            ***
*     ******                                                      ******
*     ***   Ultima revisao em: 							             *** 
*     ***   Desenvolvedor:    Tiago Tambonis                         ***
*     ***   email:            tambonis@yahoo.com.br                  ***
*     ***   msn e ymessenger: ttambonis@hotmail.com                  ***
*     ******                                                      ******
*     ***    ARQUIVOS DE ENTRADA E SAIDA MENCIONADOS NO PROGRAMA:    ***
*     ***  parametros.dat -> Paramentros de entrada                  ***
*     ***  caixamc.dat-> rede gerada com a persistência indicada     ***
*     ***                                                            ***
*     *** redesimulada.dat-> Arq de saida com a rede simulada final  ***
*     ***                                                            ***
*     ***          OS ARQUIVOS DE LEITURA E SAIDA DEVEM ESTAR NO     ***
*     ***                    DIRETORIO CORRENTE                      ***  
*     ******************************************************************
*     *  Este programa e software livre; voce pode redistribui-lo e/ou *
*     * modifica-lo sob os termos da Licenca Publica Geral GNU, con-   *
*     * forme publicada pela Free Software Foundation tanto a versao 2 *
*     * da Licenca como qualquer versao mais nova.                     *
*     *                                                                *    
*     *  Este programa e distribuido na expectativa de ser util, mas   *
*     * SEM QUALQUER GARANTIA, sem mesmo a garantia implicita de       *
*     * COMERCIALIZACAO ou de ADEQUACAO A QUALQUER PROPOSITO EM        *
*     * PARTICULAR. Consulte a Licenca Publica Geral GNU para obter    * 
*     * mais detalhes.                                                 *
*     *                                                                *  
*     * Voce deve ter recebido uma copia da Licenca Publica Geral GNU  *
*     * junto com este programa; se nao, escreva para a Free Software  *
*     * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA       *
*     * 02111-1307, USA.                                               *
*     ******************************************************************
*     *                                                              ***
*     *****************************************************************/

#include <iostream> // Fluxo de I/O por monitor, teclado...
#include <fstream>  /* i/ofstream | fstream é um manipulador de fluxos de dados de arquivos de computador
especializado para o tipo de dado nativo char. Ele permite ler e escrever em modo de texto (utiliza-se os operadores de
deslocamento de bits, << e >>) ou binário (utiliza-se os métodos read e write para buffers de dado). Fluxo de Arquivos*/
//#include <cmath>
#include <cmath>
#include <stdlib.h>
#include "MersenneTwister.h" // members of class TRandomMersenne
#include <cstdlib>

using namespace std;

//Atenção com o controle da variável bead.

class simulacaopolimero : private MTRand {
	
	private:

	   static const short int kxrede=100, kyrede=100, kzrede=100, seed=1234; //rede
	   short int xrec[100][2668],yrec[100][2668],zrec[100][2668],caixa[kxrede][kyrede][kzrede];
	   short int ncadeias,cadeia,mono,torsoes,eg,x[100], y[100], z[100], xt[100], yt[100], zt[100];
	   short int sorteiocadeia[100], maxligacoes;
	   int ncorridas,ncorridastotal;
	   float beads, A, energia_cadeia[100],temperatura,aresta, txreptation;
	   double energia_velha, energia_nova;
	   bool xestouro, yestouro, zestouro, xdirecao, ydirecao,zdirecao;
	   
	public:
	  simulacaopolimero(int seed=2697) : MTRand(seed){};
      ~simulacaopolimero(){};
	  double randon()
	  {
			return (MTRand::rand());
	   };
	  void recuperacadeias(), zeracaixa(), limpatela(), zerapolimero();
	  void zera_sorteio_cadeia(),imprime();
	  void recupera_parametros_sim(), sobreposicao(), transformacao();
	  void destransformacao(short int[], short int[], short int[]);
	  void calc_ligacoes (short int[], short int[], short int[], int[]);
	  void calc_ligacoes_totais(int [], int[]);
	  void zera_ligacoes_totais(int[]);
	  int inicializa_simulacao(),corridas(),metropolis(), vol_exluido();
	  int passacadeia(); 
      int end_move (short int, short int[], short int[], short int[], bool);
      int corner_move (short int, short int[], short int[], short int[], bool);
      int reptation_move(short int, short int[], short int[], short int[]);
      int crankshaft_move (short int, short int[], short int[], short int[],bool);
      double calc_energia_torsoes (short int[], short int[], short int[]);
      double calc_energia_vizinhos (void);
};

int main()
{	
	simulacaopolimero simulacao;
	simulacao.corridas();
	return 0;
}

/*********************************************************************
**
** Autor: Tiago Tambonis	
** Propositos: Gerenciar as corridas e obter dados apropriados.
** Modificacoes: 
** OK.
*********************************************************************/
int simulacaopolimero::corridas()
{
	
	/*O número máximo de ligacoes pi é definido no escopo de declaração 
	acima.
	* Número máximo de monomeros 100, 0 a 99.
	* 
	* Aresta máxima 0 79 -> 80. */
	
	this->recupera_parametros_sim(); //recupera os parametros da simulação
		
	this->inicializa_simulacao();
			
	return(0);
}

/*********************************************************************
**
** Autor: Tiago Tambonis.	
** Propositos: gerenciar passos adequados.
** Modificacoes: muitas... 
** 
*********************************************************************/
int simulacaopolimero::inicializa_simulacao()
{
	//beads=40, cuidado com isso.
	//ncadeias fica setado a partir da rotina recuperacadeias();
	
	static short int xtent[100], ytent[100], ztent[100];
	bool flag_end=1, flag_metropolis,flag_corner, flag_crankshaft=1;
	bool executado, flag_transf, flag_reptation;	
	short int scadeia;
	int t,cont, ligacoes[maxligacoes], ligacoestotais[maxligacoes];
	float k;
			
	MTRand mtrand1; //declaração MersenneTwister.h.
	mtrand1.seed(seed); //semente.
	
	ofstream ienergia; 
    ienergia.open("energias.dat"); 
    t=0;
    
	/* flags:
	 * 0: aceito;
	 * 1: recusado.
	*/
		
	this -> recuperacadeias(); /*recuperar as cadeias do arquivo 
	caixamc.dat geradas pelo gcadeia.cpp.*/
	imprime(); //imprime a cadeia gerada pelo gcadeia.cpp.
			
	zera_ligacoes_totais(ligacoestotais); /*zera o vetor que armazena
											o numero total de ligações
											* torsões*/
	
	for (cadeia=0;cadeia<=(ncadeias-1);cadeia++) /* cálculo da energia antiga 
													e das ligações iniciais*/
	{
		flag_transf=this->passacadeia(); //passa a cadeia setada no for.
						
		if (flag_transf!=0) this->transformacao();
		
		calc_ligacoes(x,y,z,ligacoes);
				
		calc_ligacoes_totais(ligacoes,ligacoestotais);
											
		energia_velha = energia_velha + this->calc_energia_torsoes(x,y,z);
				
		if (flag_transf!=0) this->destransformacao(x,y,z);
							
	    energia_velha = energia_velha - this->calc_energia_vizinhos();
	    		
	}
		
	cout<<"Energia inicial: "<<energia_velha<<endl;
	
	cout<<"Ligações totais iniciais: "<<endl;
	for (short int a=0;a<=maxligacoes;a++)
	{
		cout<<a<<" "<<ligacoestotais[a]<<endl;
	}
	
	cout<<"Operando..."<<endl;
	ienergia<<t<<" "<<energia_velha<<endl;
	
	getchar();
		
	cont=10000;
	cout<<endl;
			
for (ncorridas=1;ncorridas<=ncorridastotal;ncorridas++)
{	
	if (ncorridas%cont==1)
	{
		cout<<"Corrida "<<ncorridas<<endl;
		cont=cont+10000;
		imprime();
	}	
	
	zera_sorteio_cadeia();
				
	for (short int controlador_cadeia=0; controlador_cadeia<=(ncadeias-1); controlador_cadeia++)
	{
		
		start:
		
		scadeia=mtrand1.randInt(ncadeias-1);

		if (sorteiocadeia[scadeia]==1) goto start;
		else /*cadeia ainda não utilizada, é possível então realizar os 
	    movimentos nesta cadeia. */
	    {
			sorteiocadeia[scadeia]=1;
			executado=0;
			mono=mtrand1.randInt(beads-1);
			
//----------------------------------------------------------------------
//End move
			if (mono==0 || mono==(beads-1))
			{
				
				k = this->randon();
				
				if (k>=txreptation) //end move.
				{
					executado=1;
					
					cadeia=scadeia;
				
					flag_transf=this->passacadeia(); //passando a cadeia sorteada para a x,y,z da rotina end_move.
							
					if (flag_transf!=0) this->transformacao();
								
					flag_end=this->end_move(mono,xtent,ytent,ztent,flag_transf);
					
					if (!flag_end)//movimento aceito.
					{

						for (short int a=0;a<=(beads-1);++a)
						{
							xt[a]=xrec[scadeia][a];
							yt[a]=yrec[scadeia][a];
							zt[a]=zrec[scadeia][a];
							caixa[xrec[scadeia][a]][yrec[scadeia][a]][zrec[scadeia][a]]=0; //desocupo as posições que a cadeia antiga ocupava na caixa temporariamente.
						}
						for (short int a=0;a<=(beads-1);++a)
						{
							xrec[scadeia][a]=xtent[a];
							yrec[scadeia][a]=ytent[a];
							zrec[scadeia][a]=ztent[a];
							caixa[xrec[scadeia][a]][yrec[scadeia][a]][zrec[scadeia][a]]=1; //ocupo as posições que a cadeia nova ocupará na caixa.
						}	
		
						energia_nova=0.0;
						for (cadeia=0;cadeia<=(ncadeias-1);cadeia++) //cálculo da energia nova.
						{
							flag_transf=this->passacadeia();
						
							if (flag_transf!=0) this->transformacao();

							energia_nova = energia_nova + this->calc_energia_torsoes(x,y,z);
						
							if (flag_transf!=0) this -> destransformacao(x,y,z);
						
							energia_nova = energia_nova - this->calc_energia_vizinhos();		
						}
			
						flag_metropolis=this->metropolis();
						if (!flag_metropolis) //aceito a energia e a caixa permanece como está após a modificação da inserção da nova cadeia.
						{
							energia_velha=energia_nova;
							t=t+1;
							ienergia<<t<<" "<<energia_nova<<endl;
						}
						else
						{
							energia_nova=energia_velha;
							t=t+1;
							ienergia<<t<<" "<<energia_nova<<endl;
							for (short int a=0;a<=(beads-1);++a) // retornando a configuração da cadeia antiga.
							{
								caixa[xrec[scadeia][a]][yrec[scadeia][a]][zrec[scadeia][a]]=0; //desocupo as posições que a cadeia nova ocupava na caixa.
							}
							for (short int a=0;a<=(beads-1);++a) 
							{
								xrec[scadeia][a]=xt[a];
								yrec[scadeia][a]=yt[a];
								zrec[scadeia][a]=zt[a];
								caixa[xrec[scadeia][a]][yrec[scadeia][a]][zrec[scadeia][a]]=1;//ocupo as posições que a cadeia antiga ocupava na caixa.
							}
						}	
					}
				} //fim do end move
				
								
				if (k<txreptation) //reptation.
				{
					executado=1;
										
					cadeia=scadeia;
				
					flag_transf=this->passacadeia(); //passando a cadeia sorteada para a x,y,z da rotina end_move.
							
					if (flag_transf!=0) this->transformacao();
													
					flag_reptation=this->reptation_move(mono,xtent,ytent,ztent);
										
					if (!flag_reptation)//movimento aceito.
					{
						getchar();
						
						for (short int a=0;a<=(beads-1);++a)
						{
							xt[a]=xrec[scadeia][a];
							yt[a]=yrec[scadeia][a];
							zt[a]=zrec[scadeia][a];
							caixa[xrec[scadeia][a]][yrec[scadeia][a]][zrec[scadeia][a]]=0; //desocupo as posições que a cadeia antiga ocupava na caixa temporariamente.
						}
						for (short int a=0;a<=(beads-1);++a)
						{
							xrec[scadeia][a]=xtent[a];
							yrec[scadeia][a]=ytent[a];
							zrec[scadeia][a]=ztent[a];
							caixa[xrec[scadeia][a]][yrec[scadeia][a]][zrec[scadeia][a]]=1; //ocupo as posições que a cadeia nova ocupará na caixa.
						}	
		
						energia_nova=0.0;
						for (cadeia=0;cadeia<=(ncadeias-1);cadeia++) //cálculo da energia nova.
						{
							flag_transf=this->passacadeia();
						
							if (flag_transf!=0) this->transformacao();

							energia_nova = energia_nova + this->calc_energia_torsoes(x,y,z);
						
							if (flag_transf!=0) this -> destransformacao(x,y,z);
						
							energia_nova = energia_nova - this->calc_energia_vizinhos();		
						}
			
						flag_metropolis=this->metropolis();
						if (!flag_metropolis) //aceito a energia e a caixa permanece como está após a modificação da inserção da nova cadeia.
						{
							energia_velha=energia_nova;
							t=t+1;
							ienergia<<t<<" "<<energia_nova<<endl;
						}
						else
						{
							energia_nova=energia_velha;
							t=t+1;
							ienergia<<t<<" "<<energia_nova<<endl;
							for (short int a=0;a<=(beads-1);++a) // retornando a configuração da cadeia antiga.
							{
								caixa[xrec[scadeia][a]][yrec[scadeia][a]][zrec[scadeia][a]]=0; //desocupo as posições que a cadeia nova ocupava na caixa.
							}
							for (short int a=0;a<=(beads-1);++a) 
							{
								xrec[scadeia][a]=xt[a];
								yrec[scadeia][a]=yt[a];
								zrec[scadeia][a]=zt[a];
								caixa[xrec[scadeia][a]][yrec[scadeia][a]][zrec[scadeia][a]]=1;//ocupo as posições que a cadeia antiga ocupava na caixa.
							}
						}
					}
				}	
			}
//Fim do end move.
//----------------------------------------------------------------------

//----------------------------------------------------------------------	
//Corner move
			if (executado==0)
			{
				
				cadeia=scadeia;
				
				flag_transf=this->passacadeia(); /*passando a cadeia 
							sorteada para a x,y,z da rotina end_move. */
											
				if (flag_transf!=0) this->transformacao();
								
				flag_corner=this->corner_move(mono,xtent,ytent,ztent,flag_transf);
								
				if (!flag_corner) //movimento aceito
				{
								
					executado=1;
					for (short int a=0;a<=(beads-1);++a)
					{
						xt[a]=xrec[scadeia][a];
						yt[a]=yrec[scadeia][a];
						zt[a]=zrec[scadeia][a];
						caixa[xrec[scadeia][a]][yrec[scadeia][a]][zrec[scadeia][a]]=0; //desocupo as posições que a cadeia antiga ocupava na caixa.
					}
					for (short int a=0;a<=(beads-1);++a)
					{
						xrec[scadeia][a]=xtent[a];
						yrec[scadeia][a]=ytent[a];
						zrec[scadeia][a]=ztent[a];
						caixa[xrec[scadeia][a]][yrec[scadeia][a]][zrec[scadeia][a]]=1; //ocupo as posições que a cadeia nova ocupará na caixa.
					}
		
					energia_nova=0.0;
					for (cadeia=0;cadeia<=(ncadeias-1);cadeia++) //cálculo da energia nova.
					{
						flag_transf=this->passacadeia();
						
						if (flag_transf!=0) this->transformacao();

						energia_nova = energia_nova + this->calc_energia_torsoes(x,y,z);
						
						if (flag_transf!=0) this -> destransformacao(x,y,z);
						
						energia_nova = energia_nova - this->calc_energia_vizinhos();				
					}
				
					flag_metropolis=this->metropolis();
					if (!flag_metropolis) //aceito a energia e a caixa permanece como está após a modificação da inserção da nova cadeia.
					{
						energia_velha=energia_nova;
						t=t+1;
						ienergia<<t<<" "<<energia_nova<<endl;
					}
					else
					{
						energia_nova=energia_velha;
						t=t+1;
						ienergia<<t<<" "<<energia_nova<<endl;
						for (short int a=0;a<=(beads-1);++a) // retornando a configuração da cadeia antiga.
						{
							caixa[xrec[scadeia][a]][yrec[scadeia][a]][zrec[scadeia][a]]=0; //desocupo as posições que a cadeia nova ocupava na caixa.
						}
						for (short int a=0;a<=(beads-1);++a) 
						{
							xrec[scadeia][a]=xt[a];
							yrec[scadeia][a]=yt[a];
							zrec[scadeia][a]=zt[a];
							caixa[xrec[scadeia][a]][yrec[scadeia][a]][zrec[scadeia][a]]=1;//ocupo as posições que a cadeia antiga ocupava na caixa.
						}
					}	
				}
			}
//Fim do corner move.	
//----------------------------------------------------------------------

//----------------------------------------------------------------------
//Crankshat move.
			if (executado==0)
			{
				cadeia=scadeia;
				
				flag_transf=this->passacadeia(); /*passando a cadeia 
							 sorteada para a x,y,z da rotina end_move.*/
								
				if (flag_transf!=0) this->transformacao();
								
				flag_crankshaft=this->crankshaft_move(mono,xtent,ytent,ztent,flag_transf);
												
				if (!flag_crankshaft)
				{
														    														
					executado=1;
					for (short int a=0;a<=(beads-1);++a)
					{
						xt[a]=xrec[scadeia][a];
						yt[a]=yrec[scadeia][a];
						zt[a]=zrec[scadeia][a];
						caixa[xrec[scadeia][a]][yrec[scadeia][a]][zrec[scadeia][a]]=0; //desocupo as posições que a cadeia antiga ocupava na caixa.
					}
					for (short int a=0;a<=(beads-1);++a)
					{
						xrec[scadeia][a]=xtent[a];
						yrec[scadeia][a]=ytent[a];
						zrec[scadeia][a]=ztent[a];
						caixa[xrec[scadeia][a]][yrec[scadeia][a]][zrec[scadeia][a]]=1; //ocupo as posições que a cadeia nova ocupará na caixa.
					}
					
					energia_nova=0.0;
					for (cadeia=0;cadeia<=(ncadeias-1);cadeia++) //cálculo da energia nova.
					{
						flag_transf=this->passacadeia();
						
						if (flag_transf!=0) this->transformacao();

						energia_nova = energia_nova + this->calc_energia_torsoes(x,y,z);
						
						if (flag_transf!=0) this -> destransformacao(x,y,z);
						
						energia_nova = energia_nova - this->calc_energia_vizinhos();
					}
					
					flag_metropolis=this->metropolis();
					if (!flag_metropolis) //aceito a energia e a caixa permanece como está após a modificação da inserção da nova cadeia.
					{
						energia_velha=energia_nova;
						t=t+1;
						ienergia<<t<<" "<<energia_nova<<endl;
					}
					else
					{
						energia_nova=energia_velha;
						t=t+1;
						ienergia<<t<<" "<<energia_nova<<endl;
						for (short int a=0;a<=(beads-1);++a) // retornando a configuração da cadeia antiga.
						{
							caixa[xrec[scadeia][a]][yrec[scadeia][a]][zrec[scadeia][a]]=0; //desocupo as posições que a cadeia nova ocupava na caixa.
						}
						for (short int a=0;a<=(beads-1);++a) 
						{
							xrec[scadeia][a]=xt[a];
							yrec[scadeia][a]=yt[a];
							zrec[scadeia][a]=zt[a];
							caixa[xrec[scadeia][a]][yrec[scadeia][a]][zrec[scadeia][a]]=1;//ocupo as posições que a cadeia antiga ocupava na caixa.
						}
					}
				}
			}
		}//fim do for do sorteio da cadeia.
	}// fim do controlador da cadeia.
}// fim do for que controla o número de corridas.

cout<<"Energia final: "<<energia_nova<<endl;

//----------------------------------------------------------------------
//Calculo das ligações e torções.

zera_ligacoes_totais(ligacoestotais);

for (cadeia=0;cadeia<=(ncadeias-1);cadeia++) 
{
	flag_transf=this->passacadeia(); 
	
	if (flag_transf!=0) this->transformacao();
		
	calc_ligacoes(x,y,z,ligacoes);
				
	calc_ligacoes_totais(ligacoes,ligacoestotais);
}

cout<<"Ligações finais: "<<endl;
for (short int a=0;a<=maxligacoes;a++)
{
	cout<<a<<" "<<ligacoestotais[a]<<endl;
}

//----------------------------------------------------------------------
	
ienergia.close();
sobreposicao();//analisa se houve sobreposição.
imprime();//imprime a última cadeia devido ao simulacao.cpp.
cout<<"OK.";
return 0;
}// fim da função inicializa_simulacao().
/*********************************************************************
**
** Autor: Tiago Tambonis	
** Propositos: zerar as variaveis que guardam as coordenadas e as 
* posições dos monomeros das redes.
** Modificacoes: 
** OK.
*********************************************************************/
void simulacaopolimero::zerapolimero()
{
	for (short int a=0;a<=(ncadeias-1);a++)
	{	
		for (short int b=0;b<=(beads-1);b++)
		{
			xrec[a][b]=0.0;
		}
	}
}
/*********************************************************************
**
** Autor: Tiago Tambonis	
** Propositos: Limpar a tela.
** Modificacoes: 
** OK.
*********************************************************************/
void simulacaopolimero::limpatela()
{
	cout<<endl<<endl<<endl<<endl<<endl<<endl<<endl<<endl<<endl<<endl;
	cout<<endl<<endl<<endl<<endl<<endl<<endl<<endl<<endl;
}

/*********************************************************************
**
** Autor: Tiago Tambonis	
** Propositos: zerar a caixa principal.
** Modificacoes: 
** OK.
*********************************************************************/
void simulacaopolimero::zeracaixa()
{
	for (short int a=0;a<=(kxrede-1);a++)
	{
		for (short int b=0;b<=(kyrede-1);b++)
		{
			for (short int c=0;c<=(kzrede-1);c++)
			{
				caixa[a][b][c]=0.0;
			}
		}
	}
}
/*********************************************************************
**
** Autor: Tiago Tambonis	
** Propositos: Recuperar a cadeia gerada no arquivo caixamc.dat, 
* atribuindo tais valores as matrizes caixa (que vai ser considerada
* como a caixa que conterá as configurações aceitas) e a matriz caixat
* que conterá a configuração de teste no critério de metrópolis.
** Modificacoes: 
** OK.
*********************************************************************/

void simulacaopolimero::recuperacadeias()
{
	short int monomeros;
		
	ifstream caixamc; // cria o objeto cadeias.
    caixamc.open("caixamc.dat"); //abre o arquivo caixamc.dat, que deve estar no diretório corrente.
	ifstream parametros; // cria o objeto leitura.
    parametros.open("parametros.txt"); //abre o arquivo parametros.txt, que deve estar no diretório corrente.
	
//------------------------------------------------------------------------------------------------------------------------------------    
    // Leitura das propriedades do polímero.

	parametros.width(3);
    parametros>>beads; //necessidade pela compilação, não sei porque.
    parametros.width(3);
    parametros>>beads;
    parametros.close(); //fecha o arquivo setado pelo objeto desenvolvimento.
//-------------------------------------------------------------------------------------------------------------------------------------		
// Recuperação das cadeias geradas no arquivo caixamc.dat
	
	ncadeias=0;
	monomeros=0;
	
	this-> zeracaixa();
	this-> zerapolimero();
	
	while (caixamc.eof()==0.0)
	{		
		caixamc.width(3);
		caixamc>>xrec[ncadeias][monomeros];
		caixamc.width(3);
		caixamc>>yrec[ncadeias][monomeros];
		caixamc.width(3);
		caixamc>>zrec[ncadeias][monomeros];
	    caixa[xrec[ncadeias][monomeros]][yrec[ncadeias][monomeros]][zrec[ncadeias][monomeros]]=1;	

		monomeros++;
		
		if (beads==monomeros) /*Se a cadeia já foi preenchida com todos os beads uma nova deve-se iniciar.*/
		{
			ncadeias++;
			monomeros=0;
		}
	}	
	caixamc.close();
//Fim da recuperação das cadeias.
//----------------------------------------------------------------------------------------------------------------------------------
}
/*********************************************************************
**
** Autor: Ronaldo Jnio de Oliveira
** Propositos: Testar/realizar o movimento de end_move
** Modificacoes: 12/04/11, Tiago Tambonis.
** OK.
*********************************************************************/
int simulacaopolimero::end_move(short int mono, 
                                 short int xtent[], short int ytent[],
                                 short int ztent[],bool flag_transf)
{									 
   short int aleat;
   
   aleat=(int)(this->randon()*4);
	
   if(mono==0){    //�o 1o monomero
      if(x[mono]==x[mono+1] && z[mono]==z[mono+1]) // esta no eixo y
         switch (aleat){ cout<<1;
            case 0:
               xtent[mono]=x[mono+1]+1;
               ytent[mono]=y[mono+1];
               ztent[mono]=z[mono+1];
               break;
            case 1:
               xtent[mono]=x[mono+1]-1;
               ytent[mono]=y[mono+1];
               ztent[mono]=z[mono+1];
               break;
            case 2:
               xtent[mono]=x[mono+1];
               ytent[mono]=y[mono+1];
               ztent[mono]=z[mono+1]+1;
               break;
            case 3:
               xtent[mono]=x[mono+1];
               ytent[mono]=y[mono+1];
               ztent[mono]=z[mono+1]-1;
               break;
         }         
      else if(x[mono]==x[mono+1] && y[mono]==y[mono+1]) // esta no eixo z
         switch (aleat){
            case 0:
               xtent[mono]=x[mono+1]+1;
               ytent[mono]=y[mono+1];
               ztent[mono]=z[mono+1];
               break;
            case 1:
               xtent[mono]=x[mono+1]-1;
               ytent[mono]=y[mono+1];
               ztent[mono]=z[mono+1];
               break;
            case 2:
               xtent[mono]=x[mono+1];
               ytent[mono]=y[mono+1]+1;
               ztent[mono]=z[mono+1];
               break;
            case 3:
               xtent[mono]=x[mono+1];
               ytent[mono]=y[mono+1]-1;
               ztent[mono]=z[mono+1];
               break;
         }
      else if(y[mono]==y[mono+1] && z[mono]==z[mono+1]) // esta no eixo x
         switch (aleat){cout<<3;
            case 0:
               xtent[mono]=x[mono+1];
               ytent[mono]=y[mono+1]+1;
               ztent[mono]=z[mono+1];
               break;
            case 1:
               xtent[mono]=x[mono+1];
               ytent[mono]=y[mono+1]-1;
               ztent[mono]=z[mono+1];
               break;
            case 2:
               xtent[mono]=x[mono+1];
               ytent[mono]=y[mono+1];
               ztent[mono]=z[mono+1]+1;
               break;
            case 3:
               xtent[mono]=x[mono+1];
               ytent[mono]=y[mono+1];
               ztent[mono]=z[mono+1]-1;
               break;
         }
   }
   else if(mono==(beads-1)){  //�o ultimo monomero
      if(x[mono]==x[mono-1] && z[mono]==z[mono-1])   // esta no eixo y";
         switch (aleat){cout<<4;
            case 0:
               xtent[mono]=x[mono-1]+1;
               ytent[mono]=y[mono-1];
               ztent[mono]=z[mono-1];
               break;
            case 1:
               xtent[mono]=x[mono-1]-1;
               ytent[mono]=y[mono-1];
               ztent[mono]=z[mono-1];
               break;
            case 2:
               xtent[mono]=x[mono-1];
               ytent[mono]=y[mono-1];
               ztent[mono]=z[mono-1]+1;
               break;
            case 3:
               xtent[mono]=x[mono-1];
               ytent[mono]=y[mono-1];
               ztent[mono]=z[mono-1]-1;
               break;
         }         
      else if(x[mono]==x[mono-1] && y[mono]==y[mono-1])   // esta no eixo z";
         switch (aleat){cout<<5;
            case 0:
               xtent[mono]=x[mono-1]+1;
               ytent[mono]=y[mono-1];
               ztent[mono]=z[mono-1];
               break;
            case 1:
               xtent[mono]=x[mono-1]-1;
               ytent[mono]=y[mono-1];
               ztent[mono]=z[mono-1];
               break;
            case 2:
               xtent[mono]=x[mono-1];
               ytent[mono]=y[mono-1]+1;
               ztent[mono]=z[mono-1];
               break;
            case 3:
               xtent[mono]=x[mono-1];
               ytent[mono]=y[mono-1]-1;
               ztent[mono]=z[mono-1];
               break;
         }
      else if(y[mono]==y[mono-1] && z[mono]==z[mono-1])   // esta no eixo x";
         switch (aleat){cout<<6;
            case 0:
               xtent[mono]=x[mono-1];
               ytent[mono]=y[mono-1]+1;
               ztent[mono]=z[mono-1];
               break;
            case 1:
               xtent[mono]=x[mono-1];
               ytent[mono]=y[mono-1]-1;
               ztent[mono]=z[mono-1];
               break;
            case 2:
               xtent[mono]=x[mono-1];
               ytent[mono]=y[mono-1];
               ztent[mono]=z[mono-1]+1;
               break;
            case 3:
               xtent[mono]=x[mono-1];
               ytent[mono]=y[mono-1];
               ztent[mono]=z[mono-1]-1;
               break;
         }
      }
   else
     { cout << "\nTemos Problemas, end move";
		getchar();
	}
	
   if (flag_transf!=0 )
   {
		if ( xtent[mono]>(aresta-1) ) xtent[mono] = xtent[mono] - aresta;
		if ( xtent[mono]<0 ) xtent[mono] = xtent[mono] + aresta;
		if ( ytent[mono]>(aresta-1) ) ytent[mono] = ytent[mono] - aresta;
		if ( ytent[mono]<0 ) ytent[mono] = ytent[mono] + aresta;
		if ( ztent[mono]>(aresta-1) ) ztent[mono] = ztent[mono] - aresta;
		if ( ztent[mono]<0 ) ztent[mono] = ztent[mono] + aresta;  
   }

   if(!caixa[xtent[mono]][ytent[mono]][ztent[mono]]) //movimento possivel
   {   			
	   
		for (int i=0; i<=(beads-1); ++i)
		{
			if(i!=mono)
			{				
				xtent[i]=x[i];
				ytent[i]=y[i];
				ztent[i]=z[i];    
			}
        }
        
		if (flag_transf!=0 )
		{    
			for (short int a=0;a<=(beads-1);a++) if (a!=mono)
			{
				if ( xtent[a]>(aresta-1) ) xtent[a] = xtent[a] - aresta;
				if ( xtent[a]<0 ) xtent[a] = xtent[a] + aresta;
				if ( ytent[a]>(aresta-1) ) ytent[a] = ytent[a] - aresta;
				if ( ytent[a]<0 ) ytent[a] = ytent[a] + aresta;
				if ( ztent[a]>(aresta-1) ) ztent[a] = ztent[a] - aresta;
				if ( ztent[a]<0 ) ztent[a] = ztent[a] + aresta;
			}
		}
		
	return(0);
	}
	else return(1);
}	
/*********************************************************************
**
** Autor: Tiago Tambonis	
** Propositos: Calcular a energia das torsões.
** Modificacoes: 10/04/11
** OK.
*********************************************************************/
double simulacaopolimero::calc_energia_torsoes (short int xlin[], short int ylin[], 
                                  short int zlin[])
{
		double energia;

		torsoes=0; // é necessário zerar a var. "torsoes" por cadeia pois ela calculará o n.º torções por cadeia.
		
		for (short int a=0;a<=(beads-1);a++)
		{
			 if (a>=3)
			 {
				if (xlin[a]==xlin[a-1]) // eixo x
				{
					if ( (ylin[a]==ylin[a-1]) && (ylin[a]!=ylin[a-2]) ) 
					{	 
						torsoes=torsoes+1;
					}
					// posso usar else futuramente para calcular o número de ligações pi.
					if ( (zlin[a]==zlin[a-1]) && (zlin[a]!=zlin[a-2]) ) 
					{
						torsoes=torsoes+1;
					}
				}
				if (ylin[a]==ylin[a-1]) // eixo y
				{
					if ( (zlin[a]==zlin[a-1]) && (zlin[a]!=zlin[a-2]) ) 
					{
						torsoes=torsoes+1;
					}
					if ( (xlin[a]==xlin[a-1]) && xlin[a]!=xlin[a-2] )
					{
						torsoes=torsoes+1;
					}
				}
				if (zlin[a]==zlin[a-1]) // eixo z
				{
					if ( (ylin[a]==ylin[a-1]) && (ylin[a]!=ylin[a-2]) )
					{
						torsoes=torsoes+1;
					}
					if ( (xlin[a]==xlin[a-1]) && (xlin[a]!=xlin[a-2]) )
					{
						torsoes=torsoes+1;
					
					}
				}
			}
		}	
	energia=eg*torsoes;
	return (energia);
}

/*********************************************************************
**
** Autor: Tiago Tambonis	
** Propositos: Calcular a energia das torsões.
** Modificacoes: 10/04/11
** Ok.
*********************************************************************/
double simulacaopolimero::calc_energia_vizinhos()
{
	int sv;
	double energia;
	sv = 0;
	
	for (short int a=0;a<=(beads-1);a++)
	{
		sv = sv + caixa[xrec[cadeia][a]][yrec[cadeia][a]][zrec[cadeia][a]+1] + caixa[xrec[cadeia][a]][yrec[cadeia][a]][zrec[cadeia][a]-1]+
			caixa[xrec[cadeia][a]][yrec[cadeia][a]+1][zrec[cadeia][a]] + caixa[xrec[cadeia][a]][yrec[cadeia][a]-1][zrec[cadeia][a]]+
			caixa[xrec[cadeia][a]+1][yrec[cadeia][a]][zrec[cadeia][a]] + caixa[xrec[cadeia][a]-1][yrec[cadeia][a]][zrec[cadeia][a]];
	}    
	energia=A*sv;
	return (energia);
}
	
/*********************************************************************
**
** Autor: Ronaldo Jnio de Oliveira
** Propositos: Testar/realizar o movimento de "corner".
** Modificacoes: 12/04/11, Tiago Tambonis.
** 
*********************************************************************/
int simulacaopolimero::corner_move(short int mono, 
                                    short int xtent[], short int ytent[], 
                                    short int ztent[], bool flag_transf){

   short int xs, ys, zs;
   
   short int determinante = ((x[mono]*y[mono-1]*z[mono+1] + 
                              x[mono+1]*y[mono]*z[mono-1] + 
                              x[mono-1]*y[mono+1]*z[mono])-
                             (x[mono+1]*y[mono-1]*z[mono] +
                              x[mono]*y[mono+1]*z[mono-1] +
                              x[mono-1]*y[mono]*z[mono+1]));

   if (determinante)
   {
      
      xs = x[mono-1] - x[mono] + x[mono+1];
      ys = y[mono-1] - y[mono] + y[mono+1];
      zs = z[mono-1] - z[mono] + z[mono+1];
			  
	  if (flag_transf!=0)
	  {   
	    if ( xs>(aresta-1) ) xs = xs - aresta;
		if ( xs<0 ) xs = xs + aresta;
		if ( ys>(aresta-1) ) ys = ys - aresta;
		if ( ys<0 ) ys = ys + aresta;
		if ( zs>(aresta-1) ) zs = zs - aresta;
		if ( zs<0 ) zs = zs + aresta; 
	  }
	  
      if (!caixa[xs][ys][zs]) //movimento permitido
      {                        
         xtent[mono] = xs;
         ytent[mono] = ys;
         ztent[mono] = zs;
                      
         for (int i=0; i<=(beads-1); ++i) if (i!=mono)
         {
               xtent[i] = x[i];
               ytent[i] = y[i];
               ztent[i] = z[i];
         }
         if (flag_transf!=0) 
         {
			for (short int a=0;a<=(beads-1);a++) if (a!=mono)
			{
				if ( xtent[a]>(aresta-1) ) xtent[a] = xtent[a] - aresta;
				if ( xtent[a]<0 ) xtent[a] = xtent[a] + aresta;
				if ( ytent[a]>(aresta-1) ) ytent[a] = ytent[a] - aresta;
				if ( ytent[a]<0 ) ytent[a] = ytent[a] + aresta;
				if ( ztent[a]>(aresta-1) ) ztent[a] = ztent[a] - aresta;
				if ( ztent[a]<0 ) ztent[a] = ztent[a] + aresta;
			}
         }
		return (0);
      }
	  else return (1);
   }
   else 
   {
	   return (1);
	   cout<<"Problemas, corner move."<<endl;
	   getchar();
   }
}
/*********************************************************************
**
** Autor: Tiago Tambonis
** Propositos: uso do critério de metropolis.
** Modificacoes: 
** OK.
*********************************************************************/
int simulacaopolimero::metropolis()
{
	float delta,fb;
	
	delta=energia_nova - energia_velha;
		
	if (delta<=0.0)//aceito a nova energia.
	{
		return(0);
	}
	else 
	{
		fb = exp(-delta/temperatura);
		if (this->randon()<=fb)//aceito a nova energia
		{                      
			return(0);
		}
		else 
		{
			return(1); //não aceito a nova enegia.
		}
	}
}	
/*********************************************************************
**
** Autor: Ronaldo Jnio de Oliveira
** Propositos: Testar/realizar o movimento de "crankshaft".
** Modificacoes: 12/04/11, Tiago Tambonis.
** OK.
*********************************************************************/
int simulacaopolimero::crankshaft_move (short int mono, 
                                         short int xtent[],
                                         short int ytent[], 
                                         short int ztent[], bool flag_transf){
   short int xs, ys, zs;
   short int determinante;
      
//    #Observacoes ao calculo do determinante
//    Referencia: Introducao a geometria analitica, Paulo Boulos e 
//    Ivan Camargo, Makron Books, 1997 / Cap7, pag. 72.
//      
//    #Considerando v e u vetores e p uma grandeza angular
//    ||v|| ||u|| sen p = det -> u ^ v =0 <=> u//v

   if(mono!=0 || mono!=(beads-1)){
      determinante = ((x[mono]*y[mono-1]*z[mono+1] + 
                       x[mono+1]*y[mono]*z[mono-1] + 
                       x[mono-1]*y[mono+1]*z[mono])-
                      (x[mono+1]*y[mono-1]*z[mono] +
                       x[mono]*y[mono+1]*z[mono-1] +
                       x[mono-1]*y[mono]*z[mono+1]));
      
      if(determinante){
         xs = x[mono-1] - x[mono] + x[mono+1];
         ys = y[mono-1] - y[mono] + y[mono+1];
         zs = z[mono-1] - z[mono] + z[mono+1];
                 
         if(xs==x[mono+2] && ys==y[mono+2] && zs==z[mono+2]){ //�possivel 
                                                              //o crankshaft move
                                            
               //olhar em que eixo esta [mono-1] e [mono+2]
            if(abs(x[mono-1]-x[mono+2])==1){      //esta no eixo x
                                                 //exprimindo variacoes de y e z

               if(!(y[mono]-y[mono-1]))  //testar eixo z
                  if((z[mono]-z[mono-1])>0){ //z so pode diminuir
                     ztent[mono]   = z[mono]-1;
                     ztent[mono+1] = z[mono+1]-1;
                     xtent[mono]   = x[mono-1];
                     xtent[mono+1] = x[mono+2];
                     //sorteio uma direcao para y
                     if(this->randon()<0.5){
                        ytent[mono]   = y[mono]+1;
                        ytent[mono+1] = y[mono+1]+1;
                     }
                     else{
                        ytent[mono]   = y[mono]-1;
                        ytent[mono+1] = y[mono+1]-1;
                     }
                  }
                  else{    //lembrando ser impossivel variacao de z=0
                           //z so pode aumentar
                     ztent[mono]   = z[mono]+1;
                     ztent[mono+1] = z[mono+1]+1;
                     xtent[mono]   = x[mono-1];
                     xtent[mono+1] = x[mono+2];
                     //sorteio uma direcao para y
                     if(this->randon()<0.5){
                        ytent[mono]   = y[mono]+1;
                        ytent[mono+1] = y[mono+1]+1;
                     }
                     else{
                        ytent[mono]   = y[mono]-1;
                        ytent[mono+1] = y[mono+1]-1;
                     }
                  }
               else if((y[mono]-y[mono-1])<0){  //y so pode aumentar
                  ytent[mono]   = y[mono]+1;
                  ytent[mono+1] = y[mono+1]+1;
                  xtent[mono]   = x[mono-1];
                  xtent[mono+1] = x[mono+2];
                  //sorteio uma direcao para z
                  if(this->randon()<0.5){
                     ztent[mono]   = z[mono]+1;
                     ztent[mono+1] = z[mono+1]+1;
                  }
                  else{
                     ztent[mono]   = z[mono]-1;
                     ztent[mono+1] = z[mono+1]-1;
                  }
               }
               else if((y[mono]-y[mono-1]>0)){  //y so pode diminuir
                  ytent[mono]   = y[mono]-1;
                  ytent[mono+1] = y[mono+1]-1;
                  xtent[mono]   = x[mono-1];
                  xtent[mono+1] = x[mono+2];
                  //sorteio uma direcao para z
                  if(this->randon()<0.5){
                     ztent[mono]   = z[mono]+1;
                     ztent[mono+1] = z[mono+1]+1;
                  }
                  else{
                     ztent[mono]   = z[mono]-1;
                     ztent[mono+1] = z[mono+1]-1;
                  }
               }
            }
            else if(abs(y[mono-1]-y[mono+2])==1){  //esta no eixo y
                                                  //exprimindo variacoes de x e z

               if(!(x[mono]-x[mono-1]))  //testar eixo z
                  if((z[mono]-z[mono-1])>0){  //z so pode diminuir
                     ztent[mono]   = z[mono]-1;
                     ztent[mono+1] = z[mono+1]-1;
                     ytent[mono]   = y[mono-1];
                     ytent[mono+1] = y[mono+2];
                     //sorteio uma direcao para x
                     if(this->randon()<0.5){
                        xtent[mono]   = x[mono]+1;
                        xtent[mono+1] = x[mono+1]+1;
                     }
                     else{
                        xtent[mono]   = x[mono]-1;
                        xtent[mono+1] = x[mono+1]-1;
                     }
                  }
                  else{                    //z so pode aumentar
                     ztent[mono]   = z[mono]+1;
                     ztent[mono+1] = z[mono+1]+1;
                     ytent[mono]   = y[mono-1];
                     ytent[mono+1] = y[mono+2];
                     //sorteio uma direcao para x
                     if(this->randon()<0.5){
                        xtent[mono]   = x[mono]+1;
                        xtent[mono+1] = x[mono+1]+1;
                     }
                     else{
                        xtent[mono]   = x[mono]-1;
                        xtent[mono+1] = x[mono+1]-1;
                     }
                  }
               else if((x[mono]-x[mono-1])<0){   //x so pode aumentar
                  xtent[mono]   = x[mono]+1;
                  xtent[mono+1] = x[mono+1]+1;
                  ytent[mono]   = y[mono-1];
                  ytent[mono+1] = y[mono+2];
                  //sorteio uma direcao para z
                  if(this->randon()<0.5){
                     ztent[mono]   = z[mono]+1;
                     ztent[mono+1] = z[mono+1]+1;
                  }
                  else{
                     ztent[mono]   = z[mono]-1;
                     ztent[mono+1] = z[mono+1]-1;
                  }
               }
               else if((x[mono]-x[mono-1])>0){   //x so pode diminuir
                  xtent[mono]   = x[mono]-1;
                  xtent[mono+1] = x[mono+1]-1;
                  ytent[mono]   = y[mono-1];
                  ytent[mono+1] = y[mono+2];
                  //sorteio uma direcao para z
                  if(this->randon()<0.5){
                     ztent[mono]   = z[mono]+1;
                     ztent[mono+1] = z[mono+1]+1;
                  }
                  else{
                     ztent[mono]   = z[mono]-1;
                     ztent[mono+1] = z[mono+1]-1;
                  }
               }
            }
            else if (abs(z[mono-1]-z[mono+2])==1){ //esta no eixo z
                                                   //exprimindo variacoes de x e y
            
               if(!(x[mono]-x[mono-1]))  //devo testar no eixo y
                  if((y[mono]-y[mono-1])>0){   //y so pode diminuir
                     ytent[mono]   = y[mono]-1;
                     ytent[mono+1] = y[mono+1]-1;
                     ztent[mono]   = z[mono-1];
                     ztent[mono+1] = z[mono+2];
                     //sorteio uma direcao para x
                     if(this->randon()<0.5){
                        xtent[mono]   = x[mono]+1;
                        xtent[mono+1] = x[mono+1]+1;
                     }
                     else{
                        xtent[mono]   = x[mono]-1;
                        xtent[mono+1] = x[mono+1]-1;
                     }
                  }
                  else{                      //y so pode aumentar
                     ytent[mono]   = y[mono]+1;
                     ytent[mono+1] = y[mono+1]+1;
                     ztent[mono]   = z[mono-1];
                     ztent[mono+1] = z[mono+2];
                     //sorteio uma direcao para x
                     if(this->randon()<0.5){
                        xtent[mono]   = x[mono]+1;
                        xtent[mono+1] = x[mono+1]+1;
                     }
                     else{
                        xtent[mono]   = x[mono]-1;
                        xtent[mono+1] = x[mono+1]-1;
                     }
                  }
               else if((x[mono]-x[mono-1])<0){  //x so pode aumentar
                  xtent[mono]   = x[mono]+1;
                  xtent[mono+1] = x[mono+1]+1;
                  ztent[mono]   = z[mono-1];
                  ztent[mono+1] = z[mono+2];
                  //sorteio uma direcao para y
                  if(this->randon()<0.5){
                     ytent[mono]   = y[mono]+1;
                     ytent[mono+1] = y[mono+1]+1;
                  }
                  else{
                     ytent[mono]   = y[mono]-1;
                     ytent[mono+1] = y[mono+1]-1;
                  }
               }
               else if((x[mono]-x[mono-1])>0){  //x so pode diminuir
                  xtent[mono]   = x[mono]-1;
                  xtent[mono+1] = x[mono+1]-1;
                  ztent[mono]   = z[mono-1];
                  ztent[mono+1] = z[mono+2];
                  //sorteio uma direcao para y
                  if(this->randon()<0.5){
                     ytent[mono]   = y[mono]+1;
                     ytent[mono+1] = y[mono+1]+1;
                  }
                  else{
                     ytent[mono]   = y[mono]-1;
                     ytent[mono+1] = y[mono+1]-1;
                  }
               }
            }//else if
            else
               cout << "\nProblemas";
            //testando a nova ocupacao
           
           if (flag_transf!=0)
           {
			   	if ( xtent[mono]>(aresta-1) ) xtent[mono] = xtent[mono] - aresta;
				if ( xtent[mono]<0 ) xtent[mono] = xtent[mono] + aresta;
				if ( ytent[mono]>(aresta-1) ) ytent[mono] = ytent[mono] - aresta;
				if ( ytent[mono]<0 ) ytent[mono] = ytent[mono] + aresta;
				if ( ztent[mono]>(aresta-1) ) ztent[mono] = ztent[mono] - aresta;
				if ( ztent[mono]<0 ) ztent[mono] = ztent[mono] + aresta;   
				
				if ( xtent[mono+1]>(aresta-1) ) xtent[mono+1] = xtent[mono+1] - aresta;
				if ( xtent[mono+1]<0 ) xtent[mono+1] = xtent[mono+1] + aresta;
				if ( ytent[mono+1]>(aresta-1) ) ytent[mono+1] = ytent[mono+1] - aresta;
				if ( ytent[mono+1]<0 ) ytent[mono+1] = ytent[mono+1] + aresta;
				if ( ztent[mono+1]>(aresta-1) ) ztent[mono+1] = ztent[mono+1] - aresta;
				if ( ztent[mono+1]<0 ) ztent[mono+1] = ztent[mono+1] + aresta;
			}
						
            if(!(caixa[xtent[mono]][ytent[mono]][ztent[mono]]) && 
               !(caixa[xtent[mono+1]][ytent[mono+1]][ztent[mono+1]])){
               //movimento possivel
                        
               for(int i=0; i<beads; ++i) if (i!=mono && i!=(mono+1))
               {
					 xtent[i] = x[i];
                     ytent[i] = y[i];
                     ztent[i] = z[i];
               }
               
               if (flag_transf!=0 )
			   {    
					for (short int a=0;a<=(beads-1);a++) if (a!=mono and a!=(mono+1))
					{
						if ( xtent[a]>(aresta-1) ) xtent[a] = xtent[a] - aresta;
						if ( xtent[a]<0 ) xtent[a] = xtent[a] + (aresta);
						if ( ytent[a]>(aresta-1) ) ytent[a] = ytent[a] - aresta;
						if ( ytent[a]<0 ) ytent[a] = ytent[a] + (aresta);
						if ( ztent[a]>(aresta-1) ) ztent[a] = ztent[a] - aresta;
						if ( ztent[a]<0 ) ztent[a] = ztent[a] + (aresta);
					}
				}
               return(0);
            }
            else
               return(1);
         }//�possivel crankshaft
			else
				return(1);
      }//determinante==0
		else
			return(1);
   }//mono!=1 e mono!=26
	else
		return(1);
}//metodo
/*********************************************************************
**
** Autor: Tiago Tambonis	
** Propositos: simplesmente passar os valores da cadeia indicada
* pelo sorteio da cadeia e do mono para as variáveis usadas nas rotinas
* do Ronaldo.
** Modificacoes: 21/02/2011.
** OK.
*********************************************************************/
int simulacaopolimero::passacadeia()
{

	xestouro=yestouro=zestouro=0; 
		
	for (short int a=0;a<=(beads-1);a++)
	{			
		
		x[a]=xrec[cadeia][a];
		y[a]=yrec[cadeia][a];
		z[a]=zrec[cadeia][a];
		
		if (xestouro==0) /*se não houve estouro continua-se a verificar
							por estouros nas fronteiras */
		{
			if (xrec[cadeia][a]==0 && xrec[cadeia][a+1]==(aresta-1.0))
			{
				xdirecao=0;
				xestouro=1;
			}
			if (xrec[cadeia][a]==(aresta-1.0) && xrec[cadeia][a+1]==0)
			{
				xdirecao=1;
				xestouro=1;
			}
		}
		
		if (yestouro==0)
		{
			if (yrec[cadeia][a]==0 && yrec[cadeia][a+1]==(aresta-1.0))
			{
				ydirecao=0;
				yestouro=1;
			}
			if (yrec[cadeia][a]==(aresta-1.0) && yrec[cadeia][a+1]==0)
			{
				ydirecao=1;
				yestouro=1;
			}
		}
		
		if (zestouro==0)
		{
			if (zrec[cadeia][a]==0 && zrec[cadeia][a+1]==(aresta-1.0))
			{
				zdirecao=0;
				zestouro=1;
			}
			if (zrec[cadeia][a]==(aresta-1.0) && zrec[cadeia][a+1]==0)
			{
				zdirecao=1;
				zestouro=1;
			}
		}
	}
	
	if (xestouro==0 and yestouro==0 and zestouro==0) return (0);
	else return (1);
	
}
/*********************************************************************
**
** Autor: Tiago Tambonis	
** Propositos: zerar o vetor que controlará as cadeias que podem ser
* usadas no sorteio de cadeia.
** Modificacoes: 
** OK.
*********************************************************************/
void simulacaopolimero::zera_sorteio_cadeia()
{
	for (short int a=0;a<=99;a++)
	{
		sorteiocadeia[a]=0;
	}
}
/*********************************************************************
**
** Autor: Tiago Tambonis	
** Propositos: imprimir caixa principal.
** Modificacoes: 
** OK.
*********************************************************************/
void simulacaopolimero::imprime()
{
	bool flag_transf;
	
	ofstream impressao; // cria o objeto impressao.
    impressao.open("caixasimulacao.dat"); //abre o arquivo caixasimulacao, que deve estar no diretório corrente.
    
    impressao<<"X";
    impressao.width(3);
    impressao<<"Y";
    impressao.width(3);
    impressao<<"Z"<<endl;

	for(cadeia=0; cadeia<=(ncadeias-1);cadeia++)
	{
		flag_transf=this->passacadeia();
		
		if (flag_transf!=0) this->transformacao();
		
		for (short int b=0;b<=(beads-1);b++)
		{
			impressao<<x[b];
			impressao.width(3);
			impressao<<y[b];
			impressao.width(3);
			impressao<<z[b];
			impressao<<endl;
		}
		impressao<<endl<<endl<<endl;
		
	}
	impressao.close();
}
/*********************************************************************
**
** Autor: Tiago Tambonis	
** Propositos: recuperar as informações usadas na simulação.
** Modificacoes: 
** OK.
*********************************************************************/
void simulacaopolimero::recupera_parametros_sim()
{
	float const_atracao;
	
	ifstream parametrosim; // cria o objeto paramesim.
    parametrosim.open("parametrosim.txt"); //abre o arquivo parametrosim.dat, que deve estar no diretório corrente.
    
    parametrosim.width(3);
	parametrosim>>eg;
	parametrosim.width(3);
	parametrosim>>const_atracao;
	parametrosim.width(3);
	parametrosim>>temperatura;
	parametrosim.width(9);
	parametrosim>>ncorridastotal;
	parametrosim.width(3);
	parametrosim>>aresta;
	parametrosim.width(3);
	parametrosim>>maxligacoes;
	parametrosim.width(3);
	parametrosim>>txreptation;
	
	parametrosim.close();
	
	A=eg/const_atracao;
}
/*********************************************************************
**
** Autor: Tiago Tambonis	
** Propositos: verificar a existência de sobreposições.
** Modificacoes: 
** OK.
*********************************************************************/
void simulacaopolimero::sobreposicao()
{
	int posicoes[30000], c=0,d;
	short int x,y,z;

//----------------------------------------------------------------------
//Preencher o vetor que conterá todas as coordenadas.

	for (short int a=0;a<=(ncadeias-1);a++)
	{
		for (short int b=0;b<=(beads-1);b++)
		{
			posicoes[c]=xrec[a][b];
			c=c+1;
			posicoes[c]=yrec[a][b];
			c=c+1;
			posicoes[c]=zrec[a][b];
			c=c+1;		
		}
	}
//----------------------------------------------------------------------	
//Procurar sobreposições.
	
	for(c=0;c<=((beads*ncadeias)-1);c=c+3)
	{
		x=posicoes[c];
		y=posicoes[c+1];
		z=posicoes[c+2];
		for(d=c+3;d<=((beads*ncadeias)-1);d=d+3)
		{
			if (x==posicoes[d])
			{
				if (y==posicoes[d+1])
				{
					if (z==posicoes[d+2])
					{
						cout<<x<<" "<<y<<" "<<z<<endl;
						cout<<"PROBLEMA DOS GRANDES!"<<endl;
						getchar();
					}
				}
			}
		}
	}			
}
/*********************************************************************
**
** Autor: Tiago Tambonis	
** Propositos: adequar as coordenadas das cadeias que necessitam aos
* padrões exigidos pelas rotinas dos movimentos.
** Modificacoes: 
** OK.
*********************************************************************/			
void simulacaopolimero::transformacao(void)
{
	float d;
//----------------------------------------------------------------------
//Para as coordenadas x

	/* direções:
	 * 0: esquerda;
	 * 1: direita.
	*/
	
if ( xestouro==1 )
{
	x[0]=xrec[cadeia][0];
	
	for (short int a=1;a<=(beads-1);a++)
	{
		if ( abs(xrec[cadeia][a]-xrec[cadeia][a-1])!=1.0 and (xrec[cadeia][a]-xrec[cadeia][a-1])!=0)
		{
			d = (xrec[cadeia][a-1] - xrec[cadeia][a])/aresta;
			if (d<0) d=floor(d);
			if (d>0) d= ceil(d);
			x[a] = x[a-1] + d;
		}
		else
		{
			if ( (xrec[cadeia][a]-xrec[cadeia][a-1])==0 ) x[a]=x[a-1];
			else 
			{
				d = (xrec[cadeia][a] - xrec[cadeia][a-1])/aresta;
				if (d<0) d=floor(d);
				if (d>0) d= ceil(d);
				x[a] = x[a-1] + d;
			}
		}
	}
}
	
//----------------------------------------------------------------------	
// Para as coordenadas y

if ( yestouro==1 )
{
	y[0]=yrec[cadeia][0];
	
	for (short int a=1;a<=(beads-1);a++)
	{
		if ( abs(yrec[cadeia][a]-yrec[cadeia][a-1])!=1.0 and (yrec[cadeia][a]-yrec[cadeia][a-1])!=0)
		{
			d = (yrec[cadeia][a-1] - yrec[cadeia][a])/aresta;
			if (d<0) d=floor(d);
			if (d>0) d= ceil(d);
			y[a] = y[a-1] + d;
		}
		else
		{
			if ( (yrec[cadeia][a]-yrec[cadeia][a-1])==0 ) y[a]=y[a-1];
			else 
			{
				d = (yrec[cadeia][a] - yrec[cadeia][a-1])/aresta;
				if (d<0) d=floor(d);
				if (d>0) d= ceil(d);
				y[a] = y[a-1] + d;
			}
		}
	}
}

//----------------------------------------------------------------------
//Para as coordenadas z

if ( zestouro==1 )
{
	z[0]=zrec[cadeia][0];
	
	for (short int a=1;a<=(beads-1);a++)
	{
		if ( abs(zrec[cadeia][a]-zrec[cadeia][a-1])!=1.0 and (zrec[cadeia][a]-zrec[cadeia][a-1])!=0)
		{
			d = (zrec[cadeia][a-1] - zrec[cadeia][a])/aresta;
			if (d<0) d=floor(d);
			if (d>0) d= ceil(d);
			z[a] = z[a-1] + d;
		}
		else
		{
			if ( (zrec[cadeia][a]-zrec[cadeia][a-1])==0 ) z[a]=z[a-1];
			else 
			{
				d = (zrec[cadeia][a] - zrec[cadeia][a-1])/aresta;
				if (d<0) d=floor(d);
				if (d>0) d= ceil(d);
				z[a] = z[a-1] + d;
			}
		}
	}
}

//----------------------------------------------------------------------					 
}
/*********************************************************************
**
** Autor: Tiago Tambonis	
** Propositos: verificar volume excluido.
** Modificacoes: 
** OK.
*********************************************************************/
int simulacaopolimero::vol_exluido()
{
	return (0);
	
	for (short int a=0;a<=(beads-1);a++)
	{
		if (caixa[x[a]][y[a]][z[a]]!=1) 
		{		
			return (1);
			break;
		}
	}
}
/*********************************************************************
**
** Autor: Tiago Tambonis	
** Propositos: adequar as coordenadas às condições de contorno.
** Modificacoes: 
** OK.
*********************************************************************/
void simulacaopolimero::destransformacao(short int xtent[], short int ytent[],
										 short int ztent[])
{
	for (short int a=0;a<=(beads-1);a++)
	{
		if ( xtent[a]>(aresta-1) ) xtent[a] = xtent[a] - aresta;
		if ( xtent[a]<0 ) xtent[a] = xtent[a] + aresta;
		if ( ytent[a]>(aresta-1) ) ytent[a] = ytent[a] - aresta;
		if ( ytent[a]<0 ) ytent[a] = ytent[a] + aresta;
		if ( ztent[a]>(aresta-1) ) ztent[a] = ztent[a] - aresta;
		if ( ztent[a]<0 ) ztent[a] = ztent[a] + aresta;
	}
}
/*********************************************************************
**
** Autor: Tiago Tambonis	
** Propositos: Calcular a energia das torsões.
** Modificacoes: 
** 
*********************************************************************/
void simulacaopolimero::calc_ligacoes(short int xlin[], short int ylin[],
										 short int zlin[], int ligacoes[])
{
		short int nligacoes, escalar;
		
		for (short int a=0;a<=maxligacoes;a++)
		{
			ligacoes[a]=0;
		}
		
		nligacoes=0;
		for (short int a=1;a<=(beads-1);a++)
		{
			escalar = ( (x[a]-x[a-1])*(x[a]-x[a+1]) +  (y[a]-y[a-1])*(y[a]-y[a+1]) + (z[a]-z[a-1])*(z[a]-z[a+1]) );
						
			if (abs(escalar)==1) nligacoes++;
			
			if (escalar==0)
			{
				if (nligacoes!=0) ligacoes[nligacoes+1]++;
				nligacoes=0;
				ligacoes[nligacoes+1]++;
			}
		}
		
		if (nligacoes!=0) ligacoes[nligacoes+1]++;
}

/*********************************************************************
**
** Autor: Tiago Tambonis	
** Propositos: zerar o vetor que guardará o número total de ligações e
* torções.
** Modificacoes: 
** 
*********************************************************************/
void simulacaopolimero::zera_ligacoes_totais(int ligacoes_totais[])
{
	for (short int a=0;a<=maxligacoes;a++)
	{
		ligacoes_totais[a]=0;
	}
}
/*********************************************************************
**
** Autor: Tiago Tambonis	
** Propositos: Calcular o número total de ligações e torções.
** Modificacoes: 
** 
*********************************************************************/
void simulacaopolimero::calc_ligacoes_totais(int ligacoes[], int ligacoes_totais[])
{										
	for (short int a=0;a<=maxligacoes;a++)
	{
		ligacoes_totais[a]=ligacoes_totais[a]+ligacoes[a];
	}
}
/*********************************************************************
**
** Autor: Tiago Tambonis
** Propositos: Testar/realizar o movimento de "reptation".
** Modificacoes: 
** 
*********************************************************************/
int simulacaopolimero::reptation_move (short int mono, 
                                         short int xtent[],
                                         short int ytent[], 
                                         short int ztent[])
{
	bool ocupada=0;
	float r;
	
	r = this->randon();
	
	for (short int a=0;a<=(beads-1);a++)
	{
		xtent[a]=x[a];
		ytent[a]=y[a];
		ztent[a]=z[a];
	}
		
	/* Movimento */
	if (mono==0) 
	{
		
		if ( r<=1/3.0 && ztent[mono]==ztent[mono+1] && xtent[mono]==xtent[mono+1]) //espaço para crescer y.
		{

			for (short int a=0;a<=(beads-1);a++)
			{
				ytent[a]=ytent[a]-1;
			}
		}
		
		if ( (r>1/3.0 && r<=2/3.0) && ztent[mono]==ztent[mono+1] && ytent[mono]==ytent[mono+1]) //espaço para crescer x.
		{
			for (short int a=0;a<=(beads-1);a++)
			{
				xtent[a]=xtent[a]-1;
			}
		}
		
		if ( r>2/3.0 && r<=1.0 && xtent[mono]==xtent[mono+1] && ytent[mono]==ytent[mono+1]) //espaço para crescer z.
		{
			for (short int a=0;a<=(beads-1);a++)
			{
				ztent[a]=ztent[a]-1;
			}
		}
	}
	
	if (mono==(beads-1)) //último monômero.
	{
		if ( r<=1/3.0 && ztent[mono]==ztent[mono-1] && xtent[mono]==xtent[mono-1]) //espaço para crescer y.
		{
			for (short int a=0;a<=(beads-1);a++)
			{
				ytent[a]=ytent[a] + 1;
			}
		}
		
		if ( (r>1/3.0 && r<=2/3.0) && ztent[mono]==ztent[mono-1] && ytent[mono]==ytent[mono-1]) //espaço para crescer x.
		{
			for (short int a=0;a<=(beads-1);a++)
			{
				xtent[a]=xtent[a] + 1;
			}
		}
		
		if ( r>2/3.0 && r<=1.0 && xtent[mono]==xtent[mono-1] && ytent[mono]==ytent[mono-1]) //espaço para crescer z.
		{
			for (short int a=0;a<=(beads-1);a++)
			{
				ztent[a]=ztent[a] + 1;
			}
		}
	}
	
	/*Adequação às condições periódicas de contorno*/
	for (short int a=0;a<=(beads-1);a++)
	{
		if ( xtent[a]>(aresta-1) ) xtent[a] = xtent[a] - aresta;
		if ( xtent[a]<0 ) xtent[a] = xtent[a] + aresta;
		if ( ytent[a]>(aresta-1) ) ytent[a] = ytent[a] - aresta;
		if ( ytent[a]<0 ) ytent[a] = ytent[a] + aresta;
		if ( ztent[a]>(aresta-1) ) ztent[a] = ztent[a] - aresta;
		if ( ztent[a]<0 ) ztent[a] = ztent[a] + aresta;
	}
	
	for (short int a=0;a<=(beads-1);a++)
	{
		if (caixa[xtent[a]][ytent[a]][ztent[a]]==1) ocupada=1;
	}
	
	if (ocupada==1) return (1);
	else return (0);
}
