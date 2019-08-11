/* Busqueda de vecinos usando CLL:
 * 
 * El método consiste en dividir el dominio en una cantidad de Celdas N_c, cada celda tiene una regla de interaccion con las demas que determina
 * cuantas y cuales celdas tiene adyacentes. Para simplificar esta regla, se crean una serie de "Celdas Imaginarias" o "Celdas Fantasma" alrededor 
 * del dominio, esto permite que todas las celdas siempre interactuen con 8 celdas adyacentes, completándose la configuración D2Q9 para todas las
 * celdas. Estas celdas imaginarias no tienen partículas por que están fuera del dominio, sin embargo, permite detectar si por algun error numérico 
 * las partículas salen "disparadas". Si una partícula se encuentra en una Celda Fantasma, entonces algo ha salido mal en la integración temporal
 * (Verlet) o espacial (SPH) y el programa tiene que imprimir una advertencia indicando que hubo un problema. 
 */ 
 
public class Busqueda_Vecinos 
{
	//Archivo de Parámetros Generales de la Simulación
    Parametros 	Param = new Parametros();
    
    //Archivo de Funciones Auxiliares
    Utilidades 	Util  = new Utilidades();	

	//Lista de los vecinos reales de una particula i (columna), cada particula i tiene j (filas) vecinos a una distancia h_promedio = (h_i + h_j)/2.
	int[][] vecinos_reales;	
	
	//El ancho de las celdas dada por el h_i mas grande, usado como opcion para graficar las celdas.
	double ancho_de_celda;	
	
	//Lista particulas contenidas en una celda en particular y en las celdas primeras vecinas alrededor de esa celda
	// el primer indice es del nunero de celdas y el segundo es las particulas contenidas en esa celda y sus vecinas
	int[][] vecinos_potenciales; 
								
	//Asigna celda a las particulas, del primer indice: (cero , 0), el segundo indice indica la celda de la particula j
	// del segundo indice (uno, 1), el segundo indice indica en cuantas celdas mas aparece esa particula j
	int[][] celda_de_particula;	
								

	
	public  Busqueda_Vecinos(double[] h, double[][] posicion)
	{
		//Cantidad de Partículas
		int N = h.length;	
		
		//Tamaño del dominio en X
	    double xmax = Param.xmax;	
	    
	    //Tamaño del dominio en X
		double ymax = Param.ymax;	
		
		//El coeficiente K que depende de cada Kernel
		double k_kernel = new Kernel_y_Derivadas().K_Kernel();	
		
	    ancho_de_celda	 = Math.abs(Ancho_De_Celda(h, k_kernel)); 
	    
	    // Si el tamaño del ancho o largo del dominio es menor que el soporte del kernel, solo hay una celda efectiva en esa direccion, 
	    // mas dos celdas imaginarias de los costados serian tres celdas totales
		int nx = 3, ny = 3;		
								
		//Calcula el numero de celdas en el eje x
		if(xmax > ancho_de_celda)	
		{
			double paso_x = xmax/ancho_de_celda;
			double sobrante_x = (paso_x)/(int)(paso_x);
			if(sobrante_x == 1.) sobrante_x = 0.;
			nx = (int)paso_x + (int)sobrante_x + 2;
		} 

		//Calcula el numero de celdas en el eje y
		if(ymax > ancho_de_celda)	
		{
			double paso_y = ymax/ancho_de_celda;
			double sobrante_y = (paso_y)/(int)(paso_y);
			if(sobrante_y == 1.) sobrante_y = 0.;
			ny = (int)paso_y + (int)sobrante_y + 2;
		}
		
		//numero total de celdas
		int m=nx*ny;			
		
		//el numero maximo de celdas interiores
		int limite_celdas_reales_2D = m -2*(nx+ny) + 4;	
		
			//como se relacionan las celdas
		int[][] celda_veci = Regla_Interaccion_Celdas(nx, ny, limite_celdas_reales_2D);
		
		//asigna celda a las particulas
		celda_de_particula	= Celda_De_Particula(posicion, ancho_de_celda, nx,ny); 
		
		//particulas contenidas en una celda en particular y en las celdas primeras vecinas alrededor de esa celda
		vecinos_potenciales = Vecinos_Potenciales(celda_veci, celda_de_particula, m);
		
		//vecinos reales a una distancia menor o igual a h_prom. Esta es la lista que se usa para calcular las magnitudes físicas
		vecinos_reales = Vecinos_Reales(vecinos_potenciales, celda_de_particula, h, posicion, k_kernel);
		

		
	}

	// el h efectivo para calcular el tamaño de las celdas
	public double Ancho_De_Celda(double[] h, double k_kernel)	
	{

		double h_max = Util.Maximo(h);

		if(h_max <=0)
		{
			h_max = Util.Maximo(new Condiciones_Iniciales().h);
			System.out.println("Error con la distancia de suavizado que llega al metodo de vecinos");
		}
		return k_kernel*h_max;
	}

	// Acá se crean la Regla de interacción en D2Q9
	public int[][] Regla_Interaccion_Celdas(int nx, int ny, int limite_celdas_reales_2D)	
	{
		int matrices_vecinas = 8;
		int m = Math.abs(nx*ny);

		// una de estas celdas es para guardar informacion de la transformacion de celdas imaginarias a celdas reales, la otra es para guardar las particulas 
		// que se encuentren fuera del dominio a una distancia mayor a kh
		int celda_extra = 2;
		
		//celdas inmediatamente vecinas a la celda de la particula de prueba usando la regla 2D9 y organizadas en sentido de las agujas del relog, notese que 
		//la celda central (la dela particula) esta almacenada en el vector siguiente
		int[][] celda_veci = new int[m + celda_extra][matrices_vecinas + 1];
		
		//guarda la informacion de la celda de la particula
		celda_veci[m + 1 ] = new int[limite_celdas_reales_2D];	

		//Verifica que la celda no esté en una esquina
		boolean laterales=false;
		
		//Variables Auxiliares
		int k=0;
		
		//El intervalo del ciclo for ya excluye aquellas celdas que estan en la parte superior e inferior
	    for(int q=nx+1; q < m-nx-1; q++)	
	    {
	        laterales = ((q+1)%(nx)==0) || (q%(nx)==0);
	        
	        // Las celdas del centro no estan en los laterales (ni arriba o abajo ya que no se cuentan en el ciclo)
	        if(!laterales) 
	        { 
	            celda_veci[m + 1 ][k] = q;
	            celda_veci[q][0] = q-nx;
	            celda_veci[q][1] = q-nx+1;
	            celda_veci[q][2] = q+1;
	            celda_veci[q][3] = q+nx+1;
	            celda_veci[q][4] = q+nx;
	            celda_veci[q][5] = q+nx-1;
	            celda_veci[q][6] = q-1;
	            celda_veci[q][7] = q-nx-1;
	            k++; 
	        }
	    } 
	return celda_veci;
	} 


	/* Esta ultima columna extrae las celdas que se encuentran dentro del dominio computacional del metodo celda_veci, quedando un "halo" de celdas alrededor de estas 
	* celdas, en estas celdas de alrededor,  las celdas imaginarias, podrian estar las particulas fantasmas de los bordes o las particulas que se salen del dominio y
	* aun no se corrige su posicion con el metodo "Regla de Choque", estas celdas  imaginarias solo cuentan para particulas que se encuentren hasta una distancia kh del 
	* dominio computacional), el ultimo "slot" de este vector, es una celda imaginaria adicional para aquellas particulas  que estan muy lejos del dominio computacional 
	* a una distancia mayor a kh, notese que estas particulas ya no pueden interaccionar con las particulas normales del dominio debido al soporte compacto del kernel
	*/

	public int[][] Celda_De_Particula(double[][] posicion, double ancho_de_celda, int nx, int ny)
	{
		int[][] celda_de_particula = new int[2][];
		
		/* esta celda extra representa la celda de aquellas particulas alejadas del dominio computacional mas alla de kh, si la condicion CFL del tiempo funciona, ninguna 
		 * particula deberia caer en esta celda, por lo que hay que imprimir un advertencia en el programa que si alguna particula cae en esta celda
		 */
		int celda_exterior = 1;	
		
		celda_de_particula[0] = new int[posicion.length];	// indica la celda de la particula
		celda_de_particula[1] = new int[nx*ny + celda_exterior];	// cuantas celdas hay
		
		//Variables Auxiliares
		boolean xmenin, xmenex, xmayin, xmayex, ymenin, ymenex, ymayin, ymayex;
		boolean pared_izquierda, pared_derecha, pared_superior, pared_inferior;
		boolean esquina_si, esquina_ii, esquina_sd, esquina_id;
		boolean centro, paredes, esquinas;
		
		//Indican la inclinación de las paredes exteriores del dominio, siendo de forma general un paralelogramo, sin embargo, se mantiene como un cuadrado por simplificar
   		double constante_izquierda = Param.constante_izquierda;
   		double constante_derecha   = Param.constante_derecha;
    	double constante_superior  = Param.constante_superior;
    	double constante_inferior  = Param.constante_inferior;

		// Asigna celda a las particulas
		for(int i =0; i<posicion.length; i++)
		{
			xmenin = posicion[i][0] > constante_izquierda;
			xmayin = posicion[i][0] < constante_derecha;
			ymenin = posicion[i][1] > constante_superior;
			ymayin = posicion[i][1] < constante_inferior;
			centro = xmenin && xmayin && ymenin && ymayin;
			
			// Si se encuentra dentro del dominio, se le asigna una celda de manera normal
			if(centro)	
			{	
		    celda_de_particula[0][i] = (int)(posicion[i][0]/ancho_de_celda) + nx*((int)(posicion[i][1]/ancho_de_celda) + 1)  + 1;	
			celda_de_particula[1][celda_de_particula[0][i]]++;
			}
			
			// Si la particula se encuentra intensionalmente fuera del dominio (por ser particula fantasma o por haberse salido por inercia y aun no se aplica la correccion 
			// por colision) a una distancia no mayor a h del borde, se le asigna la celda real mas cercana
			
			//paredes
			if(!centro)	
			{	
				xmenex = posicion[i][0] < constante_izquierda- ancho_de_celda;
				xmayex = posicion[i][0] > constante_derecha  + ancho_de_celda;
				ymenex = posicion[i][1] < constante_superior - ancho_de_celda;
				ymayex = posicion[i][1] > constante_inferior + ancho_de_celda;
				
				pared_izquierda = !xmenin && !xmenex && ymenin && ymayin;
				pared_derecha   = !xmayin && !xmayex && ymenin && ymayin;
				pared_superior  = !ymenin && !ymenex && xmenin && xmayin;
				pared_inferior  = !ymayin && !ymayex && xmenin && xmayin;
				paredes = pared_izquierda || pared_derecha || pared_superior || pared_inferior;

				if(pared_izquierda)
				{
					celda_de_particula[0][i] = ((int)(posicion[i][1]/ancho_de_celda) + 1 )*nx + 1;
					celda_de_particula[1][celda_de_particula[0][i]]++;
				}
				if(pared_derecha)
				{
					celda_de_particula[0][i] = ((int)(posicion[i][1]/ancho_de_celda)+2)*nx - 2;
					celda_de_particula[1][celda_de_particula[0][i]]++;
				}
				if(pared_superior)
				{
					celda_de_particula[0][i] = (int)(posicion[i][0]/ancho_de_celda) + nx + 1;
					celda_de_particula[1][celda_de_particula[0][i]]++;
				}
				if(pared_inferior)
				{
					celda_de_particula[0][i] = (int)(posicion[i][0]/ancho_de_celda) + (ny-2)*nx + 1;
					celda_de_particula[1][celda_de_particula[0][i]]++;
				}

				//esquinas
				if(!paredes) 	
				{
					esquina_si = !xmenin && !xmenex && !ymenin && !ymenex;
					esquina_ii = !xmenin && !xmenex && !ymayin && !ymayex;
					esquina_sd = !xmayin && !xmayex && !ymenin && !ymenex;
					esquina_id = !xmayin && !xmayex && !ymayin && !ymayex;
					esquinas   = esquina_si || esquina_ii || esquina_sd ||esquina_id;

					if(esquina_si)
					{
						celda_de_particula[0][i] = nx + 1;
						celda_de_particula[1][celda_de_particula[0][i]]++;
					}
					if(esquina_ii)
					{
						celda_de_particula[0][i] = (ny-2)*nx + 1;
						celda_de_particula[1][celda_de_particula[0][i]]++;
					}
					if(esquina_sd)
					{
						celda_de_particula[0][i] = 2*(nx-1);
						celda_de_particula[1][celda_de_particula[0][i]]++;
					}
					if(esquina_id)
					{
						celda_de_particula[0][i] = (ny-1)*nx - 2;
						celda_de_particula[1][celda_de_particula[0][i]]++;
					}

					//celda exterior
					if(!esquinas)	
					{
						celda_de_particula[0][i] = nx*ny;
						celda_de_particula[1][celda_de_particula[0][i]]++;
					}
				}	
			}	
		}	                  
		return celda_de_particula;
	} 

	//Son aquellas particulas que estan en la misma celda de la particula de prueba o en celdas inmediatamente vecinas, usando la regla 2D9
	public int[][] Vecinos_Potenciales(int[][] celda_veci, int[][] celda_de_particula, int m)
	{
		//contador auxiliar
		int k=0;
		
		//celdas vecinas, hay 9 al ser la configuracion D2Q9, la celda central se cuenta como "0" y las demas se organizan en sentido horario
		int c0=0,c1=0,c2=0,c3=0,c4=0,c5=0,c6=0,c7=0,c8=0;
		
		//Condicionles de ayuda para evaluar que no se sobre cuenten las particulas en las celdas
		boolean par_extra_0 = false, par_extra_1 = false, par_extra_2 = false, par_extra_3 = false, par_extra_4 = false, par_extra_5 = false, par_extra_6 = false, par_extra_7 = false, par_extra_8 = false, par_efectivo = false;
		 
		//esta celda extra representa la celda de aquellas particulas alejadas del dominio computacional mas alla de kh, si la condicion CFL del tiempo funciona, ninguna particula deberia caer en esta celda, por lo que hay que 
		// imprimir un advertencia en el programa que si alguna particula cae en esta celda
		int celda_exterior = 1;	
		
		//Contador que indica la suma de las partículas en la celda de estudio y en todas las celdas adyacentes
		int suma_particulas_en_celdas =0;
		
		
		int[][] vecinos_potenciales = new int[m][0];
		
		// contador que indica cuantas particulas estan almacenadas en el area de influencia de la celda
		int[]  contador_particula_en_celda = new int[m + celda_exterior];	
		

		for(int q:celda_veci[m +1])
	    {
			suma_particulas_en_celdas =   celda_de_particula[1][celda_veci[m + 1 ][k]] + celda_de_particula[1][celda_veci[q][0]] + celda_de_particula[1][celda_veci[q][1]] + celda_de_particula[1][celda_veci[q][2]] + celda_de_particula[1][celda_veci[q][3]] + celda_de_particula[1][celda_veci[q][4]] + celda_de_particula[1][celda_veci[q][5]] + celda_de_particula[1][celda_veci[q][6]] + celda_de_particula[1][celda_veci[q][7]];
			k++;
			vecinos_potenciales[q] = new int[suma_particulas_en_celdas];
	    }

		//Se empieza a asignar los vecinos potenciales de la celda de estudio dependiendo de las adyacentes
		for(int i = 0; i<celda_de_particula[0].length; i++)
	    {
			//celda de la particula y celdas vecinas de la celda de la particula de prueba
			c0 = celda_de_particula[0][i];

			//celda exterior nunca tendra particulas validas para ser vecinas
			if(c0 != m)	
			{
				c1 = celda_veci[celda_de_particula[0][i]][0];
				c2 = celda_veci[celda_de_particula[0][i]][1];
				c3 = celda_veci[celda_de_particula[0][i]][2];
				c4 = celda_veci[celda_de_particula[0][i]][3];
				c5 = celda_veci[celda_de_particula[0][i]][4];
				c6 = celda_veci[celda_de_particula[0][i]][5];
				c7 = celda_veci[celda_de_particula[0][i]][6];
				c8 = celda_veci[celda_de_particula[0][i]][7];

				if(vecinos_potenciales[c0].length !=  0)
				{
					par_extra_0 = contador_particula_en_celda[c0] < vecinos_potenciales[c0].length;
					if(par_extra_0)
					{
					    vecinos_potenciales[c0][contador_particula_en_celda[c0]] = i;
					    contador_particula_en_celda[c0]++;
					}
					else par_efectivo = true;
				}

				if(vecinos_potenciales[c1].length !=  0)
				{
					par_extra_1 = contador_particula_en_celda[c1] < vecinos_potenciales[c1].length;
					if(par_extra_1)
					{
					    vecinos_potenciales[c1][contador_particula_en_celda[c1]] = i;
					    contador_particula_en_celda[c1]++;
					}
					else par_efectivo = true;
				}
				
				if(vecinos_potenciales[c2].length !=  0)
				{
					par_extra_2 = contador_particula_en_celda[c2] < vecinos_potenciales[c2].length;
					if(par_extra_2)
					{
					    vecinos_potenciales[c2][contador_particula_en_celda[c2]] = i;
					    contador_particula_en_celda[c2]++;
					}
					else par_efectivo = true;
				}
				
				if(vecinos_potenciales[c3].length !=  0)
				{
					par_extra_3 = contador_particula_en_celda[c3] < vecinos_potenciales[c3].length;
					if(par_extra_3)
					{
					    vecinos_potenciales[c3][contador_particula_en_celda[c3]] = i;
					    contador_particula_en_celda[c3]++;
					}
					else par_efectivo = true;
				}
				
				if(vecinos_potenciales[c4].length !=  0)
				{
					par_extra_4 = contador_particula_en_celda[c4] < vecinos_potenciales[c4].length;
					if(par_extra_4)
					{
					    vecinos_potenciales[c4][contador_particula_en_celda[c4]] = i;
					    contador_particula_en_celda[c4]++;
					}
					else par_efectivo = true;
				}
				
				if(vecinos_potenciales[c5].length !=  0)
				{
					par_extra_5 = contador_particula_en_celda[c5] < vecinos_potenciales[c5].length;
					if(par_extra_5)
					{
					    vecinos_potenciales[c5][contador_particula_en_celda[c5]] = i;
					    contador_particula_en_celda[c5]++;
					}
					else par_efectivo = true;
				}
				
				if(vecinos_potenciales[c6].length !=  0)
				{
					par_extra_6 = contador_particula_en_celda[c6] < vecinos_potenciales[c6].length;
					if(par_extra_6)
					{
					    vecinos_potenciales[c6][contador_particula_en_celda[c6]] = i;
					    contador_particula_en_celda[c6]++;
					}
					else par_efectivo = true;
				}
				
				if(vecinos_potenciales[c7].length !=  0)
				{
					par_extra_7 = contador_particula_en_celda[c7] < vecinos_potenciales[c7].length;
					if(par_extra_7)
					{
					    vecinos_potenciales[c7][contador_particula_en_celda[c7]] = i;
					    contador_particula_en_celda[c7]++;
					}
					else par_efectivo = true;
				}
				
				if(vecinos_potenciales[c8].length !=  0)
				{
					par_extra_8 = contador_particula_en_celda[c8] < vecinos_potenciales[c8].length;
					if(par_extra_8)
					{
					    vecinos_potenciales[c8][contador_particula_en_celda[c8]] = i;
					    contador_particula_en_celda[c8]++;
					}
					else par_efectivo = true;
				}

				//Advertencia que indica que determinada partícula "i" se ha salido del dominio y está en la celda extra
				if(par_efectivo)
				{
					System.out.println("Revisar el metodo de vecinos para la particula " + i);
				}

			}
		}
		return vecinos_potenciales;
	} 

	//particulas que estan a una distancia igual o menoor de kh de la particula de prueba, es una sparce matrix dado de que cada particula no tiene el mismo numero de vecinos
	public int[][] Vecinos_Reales(int [][] vecinos_potenciales, int [][] celda_de_particula, double[] h, double[][] posicion, double k_kernel)
	{
		int N = posicion.length;
		
		//contadores auxiliares
		int k = 0, u = 0;	
		
		//El h promedio entre la distancia de suavizado de dos kernels de diferentes particulas
		double h_prom = 0;
		
		//Variables auxiliares	
		double deltay=0, deltay2=0, deltax=0, deltax2=0, suma=0, radio=0, radio2=0;
		
		// Expande un poquito mas el radio del soporte para abarcar particulas que normalmente no se abarcan por problemas de presicion 
		double sensibilidad = 1.005;
		
		//Cuanta es la cantidad maxima de vecinos reales que tiene cada particula, es necesario para definir el ancho de cada columna de la matrix anterior	
		int[] vecinos_reales_maximos = new int[N];
		
		//Verifica si en verdad se encuentra dentro del rango de $kh$
		 
		boolean[] vecinos_validos 	 = new boolean[1];
		
		int[][] vecinos_reales       = new int[N][];
		
		//Este es el numero de la celda exterior
		int m 	= celda_de_particula[1].length - 1;
		int celda, cantidad_de_vecinos_potenciales, cantidad_de_vecinos_reales;

		
		
		for(int i = 0; i<N; i++)
		{
			u = 0;
			celda = celda_de_particula[0][i];
			
			//Las particulas en la celda exterior nunca seran vecinas
			if(celda != m)	
			{
				//Ve cuantas vecinos a una distancia mejor a kh tiene la particula de prueba para hacer una Sparse Matrix
				cantidad_de_vecinos_potenciales = vecinos_potenciales[celda].length;
				
				vecinos_validos = new boolean[cantidad_de_vecinos_potenciales];
				
				for(int j:vecinos_potenciales[celda])
				{
			        h_prom = (h[i]+ h[j])*0.5;
					radio = k_kernel*h_prom;
			        radio2 = radio*radio*sensibilidad;
			        deltay = posicion[i][1]-posicion[j][1];
			        deltay2 = deltay*deltay;
			        deltax = posicion[i][0]-posicion[j][0];
			        deltax2 = deltax*deltax;
			        suma = deltax2 + deltay2;
					vecinos_validos[u] = suma < radio2;	
			        if(vecinos_validos[u])
					{
						vecinos_reales_maximos[i]++;
					} 
					u++;

				}

				//Hace la Sparse Matrix y estima los vecinos reales de la particula i
				cantidad_de_vecinos_reales = vecinos_reales_maximos[i];
				vecinos_reales[i] = new int[cantidad_de_vecinos_reales];
				u = 0; 
				k = 0;
				for(int j:vecinos_potenciales[celda])
				{
			        if(vecinos_validos[u])	
			        {
			            vecinos_reales[i][k] = j;
			            k++;
			        }  
					u++;
				}

			}

			//Las particulas en la celda exterior tienen 0 vecinos, para que no interactuen entre si
			if(celda == m)	
			{
				vecinos_reales[i] = new int[0];
			}

			//Para manejar posibles errores
			if(vecinos_reales_maximos[i] == 0)	
			{
				vecinos_reales_maximos[i] = 1;
				vecinos_reales[i] = new int[vecinos_reales_maximos[i]];
				vecinos_reales[i][0] = i;
			}
		}
		
		return vecinos_reales;
	}
}
