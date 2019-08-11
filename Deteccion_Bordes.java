import java.util.*;

public class Deteccion_Bordes
{
    Parametros Param = new Parametros();
    Utilidades Util  = new Utilidades();

    

    
	public double[][] Regiones_Borde( double[] separacion)
    {
		
		int   	 	dimension 					= Param.dimension;
		int   	 	numero_condiciones_borde 	= Param.numero_condiciones_borde; 
		int[] 	 	tipo_cond_borde 		 	= Param.tipo_cond_borde;
		double[] 	tamano 						= Param.tamano;
		double[][] 	regiones_borde 				= new double[numero_condiciones_borde][2*dimension];
		
		/*solo se cuentan las regiones de los bordes de un rectangulo, sin contar fracturas
		 * en cambios tipo_cond_borde incluye todos las condiciones de borde, incluyendo las
		 * de fractura, por lo que el indice "0" se inicializa como NaN para evitar errores
		 */
		 
		regiones_borde[0][0] = Double.NaN;
		regiones_borde[0][1] = Double.NaN;
		regiones_borde[0][2] = Double.NaN;
		regiones_borde[0][3] = Double.NaN;

		/* el primer indice indica la cantidad de regiones de condiciones de borde que se van a crear, en un 
		* rectangulo normalmente serian 4, para el segundo indice el esquema es: 0 para la mitad del ancho en 
		* x del retangulo, 1 para la mitad del ancho en y del rectangulo, 2 para la ubicacion en x, 3 para la 
		* ubicacion en y
		*/
					
		//borde izquierdo
		if(numero_condiciones_borde > 1 && tipo_cond_borde[1] != 0)
		{			
			regiones_borde[1][0] =  0.6*separacion[0];
			if(tipo_cond_borde[1] == 1 || tipo_cond_borde[1] == 3)
				regiones_borde[1][1] =  0.6*tamano[1];
			if(tipo_cond_borde[1] == 2 || tipo_cond_borde[1] == 4)
				regiones_borde[1][1] =  0.6*separacion[1];
			regiones_borde[1][2] = -0.5*tamano[0];
			regiones_borde[1][3] =  0.0;
		}

		
		//borde superior
		if(numero_condiciones_borde > 2 && tipo_cond_borde[2] != 0)
		{
			if(tipo_cond_borde[2] == 1 || tipo_cond_borde[2] == 3)
				regiones_borde[2][0] = 0.6*tamano[0];
			if(tipo_cond_borde[2] == 2 || tipo_cond_borde[2] == 4)
				regiones_borde[2][0] = 0.6*separacion[0];
			regiones_borde[2][1] = 0.5*separacion[1];
			regiones_borde[2][2] = 0.0;
			regiones_borde[2][3] = 0.5*(1.0*separacion[1] - tamano[1]);
		}
		
		//borde derecho
		if(numero_condiciones_borde > 3 && tipo_cond_borde[3] != 0)
		{
			regiones_borde[3][0] =  0.6*separacion[0];
			if(tipo_cond_borde[3] == 1 || tipo_cond_borde[3] == 3)
				regiones_borde[3][1] =  0.6*tamano[1];
			if(tipo_cond_borde[3] == 2 || tipo_cond_borde[3] == 4)
				regiones_borde[3][1] =  0.6*separacion[1];
			regiones_borde[3][2] =  0.5*tamano[0];
			regiones_borde[3][3] =  0.0;
		}	
		
		//System.out.println(tipo_cond_borde[2] + " " + regiones_borde[2][1] + " " + regiones_borde[2][3] + " " + 0.6*tamano[1] + " " + 0.6*separacion[1]);
		
		//borde inferior
		if(numero_condiciones_borde > 4 && tipo_cond_borde[4] != 0)
		{
			if(tipo_cond_borde[4] == 1 || tipo_cond_borde[4] == 3)
				regiones_borde[4][0] = 0.5*tamano[0];
			if(tipo_cond_borde[4] == 2 || tipo_cond_borde[4] == 4)
				regiones_borde[4][0] = 0.5*separacion[0];
			regiones_borde[4][1] = 0.5*separacion[1];
			regiones_borde[4][2] = 0.0;
			regiones_borde[4][3] = 0.5*(tamano[1]  - 1.0*separacion[1]);
		}
		
		//borde izquierdo, segunda condicion de borde
		if(numero_condiciones_borde > 5 && tipo_cond_borde[5] != 0)
		{			
			regiones_borde[5][0] =  0.6*separacion[0];
			if(tipo_cond_borde[5] == 1 || tipo_cond_borde[5] == 3)
				regiones_borde[5][1] =  0.6*tamano[1];
			if(tipo_cond_borde[5] == 2 || tipo_cond_borde[5] == 4)
				regiones_borde[5][1] =  0.6*separacion[1];
			regiones_borde[5][2] = -0.5*tamano[0];
			regiones_borde[5][3] =  0.0;
		}
		
		return regiones_borde;
    }
    
    public boolean[] Esta_En_borde(int particula, double[] posicion, double[][] regiones_borde)
    {
		boolean[] esta_en_borde = new boolean[regiones_borde.length];
		if(!(posicion[0] == 0 && posicion[1] ==0))
		{
			for(int pared = 1; pared < regiones_borde.length; pared++)
			{
			esta_en_borde[pared] = (posicion[0] >= regiones_borde[pared][2] - regiones_borde[pared][0] && posicion[0] <= regiones_borde[pared][2] + regiones_borde[pared][0]) && (posicion[1] >= regiones_borde[pared][3] - regiones_borde[pared][1] && posicion[1] <= regiones_borde[pared][3] + regiones_borde[pared][1]);
			//if(particula == 71 && pared == 2) System.out.println(esta_en_borde[pared] + " " + posicion[1] + " " + (regiones_borde[pared][3] - regiones_borde[pared][1]) + " " + (regiones_borde[pared][3] + regiones_borde[pared][1]));
			//if(particula == 71 && pared == 2) System.out.println(esta_en_borde[pared] + " " + posicion[1] + " " + regiones_borde[pared][3] + " " + regiones_borde[pared][1]);
			//System.out.println(particula + " " + pared + " " + posicion[0] + " " + (regiones_borde[pared][2] - regiones_borde[pared][0]) + " " + (regiones_borde[pared][2] + regiones_borde[pared][0]));
			}
		
		}
		return esta_en_borde;
	}
	

	
	public int[] borde_de_particula(List<Integer> particulas_fractura, List<Integer> particulas_en_borde, double[] separacion, double[][] posicion)
	{
		/* Este metodo clasifica las particulas si estan en bordes de un rectangulo, regresa -2 si la particula esta en el centro, -1 si esta en una
		 * fractura o cualquier otro borde, 0 si esta en la pared izquierda, 2 si esta en la derecha, 1 en la superior y 3 en la inferior
		 * */
		 
		boolean[] 	esta_en_borde 				= new boolean[1];
		boolean		sin_pared					= true;
		int 		N 							= posicion.length;
		int 		numero_condiciones_borde 	= Param.numero_condiciones_borde;
		int[] 		borde_de_particula 			= Util.Vector_Componentes_Iguales(N, -2);
		int[] 		orden_condiciones_borde 	= Param.orden_condiciones_borde;
		double[][] 	regiones_borde 				= Regiones_Borde(separacion);
		int nx				= Param.nx;
		int ny				= Param.ny;
		int n_sobrante 		= Param.n_sobrante;
		int a = 1;
		
		
		

		
		//Util.Imprimir_Lista(particulas_fractura);
		
		/* se selecionan los bordesde las particulas dependiendo de que pared esten, notese que si una pared está en 
		* un borde y a la vez en una fractura, se tiene prioridad la clasificacion del borde_de_particula
		*/
		for(int i : particulas_en_borde)
		{
			sin_pared 	  = true;
			esta_en_borde = Esta_En_borde(i, posicion[i], regiones_borde);
			
			//System.out.println(i);
			
			
			for(int pared: orden_condiciones_borde)	
			{
				/* los indices de esta_en_borde comienza en 0 para la pared izquierda y termina en 3 para la inferior,
				 *  en cambio, los indices de orden_condiciones_borde y borde_de_particula tiene que ser 0 para
				 * fracturas, 1 para izquierda, 4 para inferior, es decir, estan añadidos uno mas, por lo tanto hay
				 * que restarle uno al indice "pared" para que comience igual que los indices de esta_en_borde
				 */
				 
			//System.out.println(i + " " + pared);
				if(pared > 0 && esta_en_borde[pared]) 
				{
					//if(i== Param.particula_prueba) 
					
					 borde_de_particula[i] = pared;
					 sin_pared = false;
					
				}
				if(pared == orden_condiciones_borde[orden_condiciones_borde.length - 1] && sin_pared)
					borde_de_particula[i] = -1;	// Si no se le encuentra condicion de borde se le coloca en el indice -1
			}

		}
		
		//incluye las particulas que esten en particulas_fractura como particulas de los bordes
		for(int i : particulas_fractura)
		{
			borde_de_particula[i] = 0;	
		}
		
		for(int i = 0; i < borde_de_particula.length; i++)
			if(borde_de_particula[i] == -1) System.out.println(i);
		
		return borde_de_particula;
	}
}

/* Esta clase retorna la superficie de aquellas particulas que estan en el borde,
 * solo toma en cuenta las de los bordes fisicos/fracturas , por lo que debe ser 
 * corrido con el condicional de borde_de_particula[i] != -2, ya que cuando 
 * borde_de_particula[i] == -2 entonces esa  particula "i" es del bulk y no tiene
 * superficie fisica del medio
 */
 
class Superficie_Particulas
{
	/* este array tiene que tener la longitud de particulas_fractura + 
	 * particulas_en_borde menos aquellas particulas que esten en ambas listas
	 * y hay que sacar una de sus copias, de todas formas esto es igual a la
	 * cantidad de veces que ocurre la condicion borde_de_particula[i] != -2
	*/
	double[] superficie;	
	Utilidades Util = new Utilidades();
    Parametros Param = new Parametros();   

	int			N = Param.N;
	int 		dimension 	= Param.dimension;
			    
	public Superficie_Particulas(List<Integer> particulas_en_borde, int[][] vecinos, double[][] posicion, boolean[][] vecinos_hook)
	{
	
		double		distancia	= 0;
		superficie 	= new double[N];
		double angulo;
		double angulo_tolerancia = 45./180.*Param.PI;
		double distancia_media = Param.separacion;
		double[] vector_1 = new double[dimension];
		double[] vector_2 = new double[dimension];
		int[][] vecinas_borde = Vecinas_Borde(particulas_en_borde, vecinos, posicion, vecinos_hook);
		
		for(int i: particulas_en_borde)
		{
			if(vecinas_borde[i][0] != -1 && vecinas_borde[i][1] != -1)
			{
				vector_1 = Util.Resta_Vectores(posicion[i], posicion[vecinas_borde[i][0]]);
				vector_2 = Util.Resta_Vectores(posicion[i], posicion[vecinas_borde[i][1]]);
				angulo = Util.Angulo_2D(vector_1, vector_2);
				distancia_media = 0.5*(Util.Modulo_Vector(vector_1) + Util.Modulo_Vector(vector_2));
				
				if(Math.abs(angulo) > angulo_tolerancia)
					superficie[i] = Param.RAIZ_DOS*distancia_media;
				else
					superficie[i] = distancia_media;
			}
			else
				superficie[i] = distancia_media;
		}
		//Util.Imprimir_Lista(particulas_en_borde);
		//Util.Imprimir_Vector(superficie);
				
	}
	
	//este metodo busca las dos particulas del borde mas cercanas a una particula en el borde
	public int[][] Vecinas_Borde(List<Integer> particulas_en_borde, int[][] vecinos, double[][] posicion, boolean[][] vecinos_hook)
	{
		int[][] vecinas_borde = Util.Matriz_Componentes_Iguales(N, dimension, -1);
		
		int 		vecina;
		List<Double> 	primera_distancia_mas_corta = new ArrayList<Double>();
		List<Integer> 	indice_distancia_mas_corta = new ArrayList<Integer>();
		
		//se hace un for buscando en todas las particulas que estan en los bordes (solo estas tienen
		// superficie, las demas no tienen).
		for(int i: particulas_en_borde)
		{

			//Se busca en los vecinos de estas particulas
			for(int j = 0 ; j < vecinos[i].length; j++)
			{
				
				vecina = vecinos[i][j];
				
				//pero solo los vecinos que esten en la superficie, que no sean ellas mismas y que
				//esten conectadas cohesivamente (que formen una superficie entre ellas y no exista
				//una fractura en el medio)
				if (i != vecina && particulas_en_borde.contains(vecina) && vecinos_hook[i][j])
				{
					primera_distancia_mas_corta.add(Util.Distancia(posicion[i], posicion[vecina]));
					indice_distancia_mas_corta.add(vecina);
				}
			}

			if(primera_distancia_mas_corta.size() == 0)
				System.out.println(" En el metodo Vecinas_Borde, clase superficie_fractura y archivo deteccion_bordes, hay un problema con la particula " + i + " que no tiene mas vecinos a su alrededor que pertenezcan al borde/fractura. Su cantidad de vecinos total es: " + vecinos[i].length);
			else
			{
			//se guardan los dos vecinos del borde mas cercanos a la particula i del borde
			int[] wolverine = Util.Indice_Minimo_Dos_Primeros(primera_distancia_mas_corta);
			vecinas_borde[i] = new int[]{indice_distancia_mas_corta.get(wolverine[0]), indice_distancia_mas_corta.get(wolverine[1])};
			}
			//se reinician los elementos para la siguiente iteracion
			primera_distancia_mas_corta = new ArrayList<Double>();
			indice_distancia_mas_corta = new ArrayList<Integer>();

		}
		
	return vecinas_borde;	
	}
}



class Particulas_En_Borde_obsoleto
{
	
	boolean[]  	esta_en_borde;
	int[] 		particulas_en_borde;
	double[][] 	vector_normal;
	double[][]	vector_tangente;
	double[][]	direccion_efectiva;
	
	Utilidades Util = new Utilidades();
	
	public Particulas_En_Borde_obsoleto(Parametros Param, Utilidades Util, int N, int[][] vecinos, double[] tamano, double[] separacion, double[][] posicion)
	{
		
		int			dimension	= Param.dimension;
		int 		largo 		= 0;
		int 		k 			= 0;
		int 		particula;
		int			otra_particula;
		int[]		indices_angulos;
		int[][]     solo_vecinos_distancia = Solo_Vecinos_Distancia(vecinos, 1.9*separacion[0], posicion);	
		double 		tolerancia 	= 1.5;
		double 		radian 		= 57.2857;
		double		RAIZ_DOS	= Param.RAIZ_DOS;
		double 		dx, dy, modulo;
		double 		dif_maxima	= radian;
		double[] 	diferencia 	= new double[largo];
		double[] 	angulo_sin_ordenar 	= new double[largo];
		double[] 	angulos 	= new double[largo];
		double[] 	dominio  	= Util.Vector_Componentes_Iguales(N, tolerancia*separacion[0]);
		
		
		direccion_efectiva 	= new double[N][posicion[0].length];	
		esta_en_borde 		= new boolean[N];
		vector_normal  		= new double[N][dimension];
		vector_tangente 	= new double[N][dimension];

		List<Integer> array_particulas_en_borde = new ArrayList<Integer>();

		
		if(N>1)
		{
			for(int i =0; i < N; i++)
			{
				k 		= 0;
				largo 	= solo_vecinos_distancia[i].length;
							
		
							
				angulo_sin_ordenar 	= new double[largo];
				diferencia  		= new double[largo];
				angulos 			= new double[largo];
				indices_angulos		= new int[largo];
				
				for(int j : solo_vecinos_distancia[i])
				{
					dx = posicion[j][0] - posicion[i][0];
					dy = posicion[j][1] - posicion[i][1];
					direccion_efectiva[i][0] += dx;
					direccion_efectiva[i][1] += dy;
					angulo_sin_ordenar[k] = Util.Angulo_2D(dx, dy);					
					k++;


				}
							//	Util.Imprimir_Vector(solo_vecinos_distancia[i]);
		
					
				Ordenar_Vector_Menor_Mayor Sort =  new Ordenar_Vector_Menor_Mayor(angulo_sin_ordenar);
				
				angulos 		= Sort.array_double;
				indices_angulos = Sort.indice;
					
				diferencia[0] = radian*(angulos[0] - angulos[largo - 1]);
				
				/*
				if(esta_particula_fracturada)
				{
					
					esta_en_borde[i] = true;
					array_particulas_en_borde.add(i);
					vector_normal[i][0] = vector_normal_fractura[indice_particula_fracturada][0];
					vector_normal[i][1] = vector_normal_fractura[indice_particula_fracturada][0];
					
				}
				
				else
				* */
				{
					for(int q = 0; q < largo; q++)
					{
						
						if(q > 0)
							diferencia[q] = radian*(angulos[q] - angulos[q-1]); 
						
						if(diferencia[q] >= 360)	
							diferencia[q] -= 360;
							
						if(diferencia[q] < 0)	
							diferencia[q] += 360;
						
						//System.out.println(i + " " + q + " " + diferencia[q] + " " + radian);
						
						if( diferencia[q] > radian)
						{
							esta_en_borde[i] = true;
							array_particulas_en_borde.add(i);

							if(q == 0)
							{
								particula 	   = solo_vecinos_distancia[i][indices_angulos[0]];
								otra_particula = solo_vecinos_distancia[i][indices_angulos[largo - 1]];						
							}
							else
							{
								particula 	   = solo_vecinos_distancia[i][indices_angulos[q]];
								otra_particula = solo_vecinos_distancia[i][indices_angulos[q-1]];
							}
							
							dx = posicion[particula][0] - posicion[otra_particula][0];
							dy = posicion[particula][1] - posicion[otra_particula][1];
							modulo = Math.sqrt(dx*dx + dy*dy);
							
							vector_tangente[i][0] = dx/modulo;
							vector_tangente[i][1] = dy/modulo;
							
							
							//if(i==6) //System.out.println(otra_particula + " " + vector_tangente[i][0] + " " + vector_tangente[i][1]);
							//Util.Imprimir_Vector(solo_vecinos_distancia[i]);


							vector_normal[i] = Util.Vector_Perpendicular_Horario(vector_tangente[i]);
							
							if(Util.Producto_Punto(vector_normal[i], direccion_efectiva[i]) > 0)
								vector_normal[i] = Util. Vector_Perpendicular_AntiHorario(vector_tangente[i]);
								
						}	
					}
				}		
			}
		}
		else
		{
			array_particulas_en_borde.add(0);
		}
		
		//Util.Imprimir_Matriz(vector_normal);
		particulas_en_borde = array_particulas_en_borde.stream().filter(ñ -> ñ != null).mapToInt(ñ -> ñ).toArray();

	}
	
	    
    public int[][] Solo_Vecinos_Distancia(int[][] vecinos, double separacion, double[][] posicion)
	{
		int largo = vecinos.length;
		int ancho;
		int k = 0;
		int[][] nueva_lista = new int[largo][];
		
		double distancia = 0;
		
		List<Integer> lista_auxiliar;
		
		for(int i = 0; i < largo; i++)
		{
			lista_auxiliar = new ArrayList<Integer>();
			

			for(int j: vecinos[i])
			{
				distancia = Util.Distancia(posicion[j], posicion[i]);
				
				if(i != j && distancia < separacion)
				{
					//System.out.println(i + " " + j + " " + distancia);
					lista_auxiliar.add(j);
				}
			}
			nueva_lista[i] = lista_auxiliar.stream().filter(ñ -> ñ != null).mapToInt(ñ -> ñ).toArray();
			//System.out.println(i + " " + nueva_lista[i].length + " " + distancia);
		}
		
		return nueva_lista;
	}
}
