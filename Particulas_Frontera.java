import java.util.*;
import java.util.stream.*;

public class Particulas_Frontera{}

class Particulas_En_Borde
{
	Utilidades Util = new Utilidades();
	Parametros Param = new Parametros();
	
	boolean[][] 	vecinos_hook;
	List<Integer>	particulas_en_borde;
	double[] 		superficie;
	double[][] 		vector_normal;
	double[][] 		vector_tangente;
			
		
	public Particulas_En_Borde
	(
		double separacion, 
		double[] h, 
		int[][] vecinos, 
		double[][] kernel, 
		double[][][] grad_kernel, 
		List<Integer> particulas_fractura, 
		double[] punto_origen_fractura, 
		double[] punto_borde_fractura,
		double[][] posicion,
		double masas,
		double[] densidad
	)
	{
		
		boolean cond_a;
		boolean cond_b;
		
		int N = Param.N;
		int dimension = Param.dimension;
		int vecino;
		
		double	profundidad					= 0;
		double 	RAIZ_DOS					= Param.RAIZ_DOS;
		double 	modulo 						= 0.0;
		double 	UNO_ENTRE_modulo 			= 0.0;
		double 	tolerancia_kernel 			= 1.0 - 0.1;
		double 	tolerancia_modulo 			= Math.pow(10, -5);			
		double[] kernel_distancia_raiz_dos 	= Vector_Kernel_Constante(h, Util.Vector_Componentes_Iguales(N, RAIZ_DOS*separacion), dimension);
		double[] vector_distancia_ij = new double[dimension];
		
		//Util.Imprimir_Lista(particulas_fractura);
		vecinos_hook  		= new boolean[N][1];
		//particulas_en_borde	= new ArrayList<Integer>();
		particulas_en_borde	= particulas_fractura.stream().collect(Collectors.toList());
		superficie			= new double[N];
		vector_normal 		= new double[N][dimension];
		vector_tangente 	= new double[N][dimension];
		
		for(int i = 0; i < N; i++)
		{
			vecinos_hook[i] = new boolean[vecinos[i].length];
			
			for(int j = 0; j < vecinos[i].length; j++)
			{
				vecinos_hook[i][j] = true;
				
				vecino = vecinos[i][j];
				Interseccion_Lineas Interseccion = new Interseccion_Lineas(posicion[i], posicion[vecino], punto_origen_fractura, punto_borde_fractura);
				
				if(Interseccion.se_intersectan)
					vecinos_hook[i][j] = false;
					
				cond_a = (kernel[i][j] - kernel_distancia_raiz_dos[i]) > 0;
				
				if(cond_a && vecinos_hook[i][j])
				{
					//if(grad_kernel[i][j][0] == grad_kernel[i][j][0])
					vector_normal[i][0] += -grad_kernel[i][j][0];
					
					//if(grad_kernel[i][j][1] == grad_kernel[i][j][1])
					vector_normal[i][1] += -grad_kernel[i][j][1];
					
					
				}
				//System.out.println(i + " " + vecino +  " "  + cond_a + " " + vecinos_hook[i][j] + " " + grad_kernel[i][j][0]);
			}
				
			
			modulo = Util.Modulo_Vector(vector_normal[i]);
			
			//if(i == 293)System.out.println(modulo);
			
			
			if(modulo > tolerancia_modulo && !particulas_en_borde.contains(i))
			{
				//System.out.println(i + " " + modulo);
				UNO_ENTRE_modulo = 1.0/modulo;
				particulas_en_borde.add(i);
				
				/* este ciclo es para saber cual es la particula del bulk que está indenmediatamente "debajo" de una particula del borde, esto es para calcular 
				* la "profundidad" del area de influencia de la particula del borde y se despeja la superficie a partir del "volumen" de la particula de SPH
				*  evidentemente solo es algo para 2D y distribuciones cercanas a rectangulares */
				
				for(int j = 0; j < vecinos[i].length; j++)
				{
					vecino = vecinos[i][j];
					
					vector_distancia_ij = Util.Resta_Vectores(posicion[i], posicion[vecino]);
					cond_a 				= (kernel[i][j] - kernel_distancia_raiz_dos[i]*tolerancia_kernel) > 0.0;
					cond_b 				= Util.Modulo_Vector(Util.Producto_Cruz(vector_distancia_ij, vector_normal[i])) <= tolerancia_modulo ;
				
					if(i != vecino && (cond_a && cond_b && vecinos_hook[i][j] || particulas_fractura.contains(i)))
					{
						profundidad = Util.Modulo_Vector(vector_distancia_ij);
						superficie[i] = (masas/densidad[i])/profundidad;	
						//System.out.println(i + " " + vecino + " " + superficie[i] + " " + separacion);
					}
				}
			}
			
			if(modulo < tolerancia_modulo)
				UNO_ENTRE_modulo = 0.0;
				
			vector_normal[i][0] *= UNO_ENTRE_modulo;
			vector_normal[i][1] *= UNO_ENTRE_modulo;
		
			vector_tangente[i] = Util.Vector_Perpendicular_AntiHorario(vector_normal[i]);
			
		}
		
		//System.out.println(vector_normal[168][0] + " " + vector_normal[168][1]);
		//Util.Imprimir_Lista(particulas_fractura);
	}
	
	public double Kernel_Constante(int i, double h, double separacion, int dimension)
    {
        double q = separacion/(h);
        return  new Kernel_y_Derivadas().Kernel(-1, -1, q, h, dimension);
    } 
       
	public double[] Vector_Kernel_Constante(double[] h, double[] separacion, int dimension)
	{

		double[] kernel_distancia = new double[h.length];
		
		for(int i =0; i < h.length; i++)
			kernel_distancia[i] = Kernel_Constante(i, h[i], separacion[0], dimension);
			
		return kernel_distancia;
	}
}


/* Clase que calcula si dos lineas se intercectan y cual es el punto de interseccion en caso de que las lineas SI se intercecten
 */
class Interseccion_Lineas
{
	boolean  se_intersectan;	//retorna "true" si las lineas se intersectan y "false" en caso de que no
	double[] posicion_puntos;	/* en caso de ser "true" la variable anterior, retorna las coordenadas del punto de intersección
								* AVISO: se debe hacer una comprobación de que las lineas si se intersectan antes de usar esta
								* variable, sino, retorna el punto [Nan,Nan] de manera intencionada
								*/
	
    Utilidades 	Util  = new Utilidades();
    Parametros 	Param = new Parametros();
    
		//Se intenta buscar si las linea 1 (formada por el punto 1 y 2) se intersecta con la linea 2 (formada por los puntos 3 y 4)
	public Interseccion_Lineas( double[] punto_1, double[] punto_2, double[] punto_3, double[] punto_4)
	{
		double rxs;
		double qpxr;
		double qpxs;
		double a;
		double b;
		double c;
		double tolerancia 	= -1*Math.pow(10,-5);	//para evitar problemas de redondeo y coma flotante
		
		double[] p 	=  Util.Copiar_Vector(punto_1);
		double[] q 	=  Util.Copiar_Vector(punto_3);
		double[] r 	= Util.Resta_Vectores(punto_2, punto_1);
		double[] s 	= Util.Resta_Vectores(punto_4, punto_3);
		double[] qp = Util.Resta_Vectores(q, p);
			
		rxs 	= Util.Producto_Cruz(r, s)[0];
		qpxr 	= Util.Producto_Cruz(qp, r)[0];
		qpxs 	= Util.Producto_Cruz(qp, s)[0];
		
		se_intersectan 	= false;
		posicion_puntos = Util.Vector_Componentes_Iguales(punto_1.length, Double.NaN);
		/* se inicializa el punto de interseccion con Nan para evitar el caso donde las lineas no intersectan y se regrese el valor [0,0] 
		 * sin ningun warning o error y el programa haga las cuentas como sí si existiera una intersepcion entre las lineas en el origen
		 */
			
		
		//algoritmo simplificado para 2d de  Ronald Goldman, published in Graphics Gems
		if(Math.abs(rxs) <= tolerancia)								//son paralelas
		{
			if(Math.abs(qpxr) <= tolerancia)							//son colineales, hay que verificar si una esta contenida en la otra
			{
					
				double r2 	= Util.Producto_Punto(r, r);
				
				a =  Util.Producto_Punto(qp, r)/r2;
				c =  a + Util.Producto_Punto(s, r)/r2;			
	
				b = Math.max(a, c);					//esto evita el problema de saber si s*r es mayor o menor que cero para organizar el intervalo
				a = Math.min(a, c);					//"a" es el valor menor y "b" el mayor del intervalo
				
				
				
				//todo el problema quedo reducido a interseccion de dos intervalos en 1D, el intervalo [a,b] y [0,1]
				if(a <= 0.0 || (b - 1.0) >= 0.0)
				{					
					se_intersectan = false;			//si el valor maximo es menor que cero o el minimo es mayor que uno, los intervalos no se intersectan, notese que incluye los puntos del intervalo
				
				}
				else
				{
					se_intersectan = true;			//del resto si se intersectan
					
				}
				
				//se_intersectan = false;
			}
			
			else									//son paralelas y nunca se cruzan
				se_intersectan = false;
		}
		
		else										//no son paralelas
		{
			a = qpxs/rxs;
			b = qpxr/rxs;
			
			
			
			if(a >= tolerancia && (a - 1.0) <= tolerancia && b >= tolerancia && (b - 1.0) <= tolerancia)	
			{										//no son paralelas y se cruzan, normalmente las desigualdades serian >= o <=. pero quiero que los puntos extremos No esten incluidos
				posicion_puntos[0] = p[0] + a*r[0];
				posicion_puntos[1] = q[1] + b*s[1];
				se_intersectan = true;
				//System.out.println(q+ " " + p + " " + a + " " + b + " " + se_intersectan);	
			}
			else 
			{									//no son paralelas pero no se intersectan
				se_intersectan = false;
				//System.out.println(q+ " " + p + " " + a + " " + b + " " + se_intersectan);	
			}
		}
	}
}


class Particulas_Fractura
{
	Parametros 		Param 		= new Parametros();
	Utilidades 		Util 		= new Utilidades();
	    
	double[] 		punto_origen_fractura;
	double[]		punto_borde_fractura;
	double[]		punto_origen_fractura_recortado;
	double[]		extremos_rectangulo;
	double[][] 		rectangulo_particulas;
	List<Integer> 	particulas_fractura;
	
	public Particulas_Fractura(double[] centro_factura_recta, double[] tamano_fractura, double[][] posicion_particulas, double separacion)
	{
		double PI 						= Param.PI;
		double ancho_rectangulo			= 1.1*separacion;
		double extra_largo_rectangulo 	= 0.0*separacion;
		double angulo_fractura 			=  tamano_fractura[1]/180.0*PI;
		double tolerancia_longitud		= 0.1;	//para evitar problemas de presicion al borde de las fractura
		double longitud_fractura_ini	= tamano_fractura[0] + tolerancia_longitud*separacion;
			
		punto_borde_fractura 	= Punto_Borde_Fractura(centro_factura_recta, longitud_fractura_ini, angulo_fractura);
		punto_origen_fractura 	= Punto_Origen_Fractura(centro_factura_recta, longitud_fractura_ini, angulo_fractura);
		rectangulo_particulas	= Rectangulo(punto_origen_fractura, punto_borde_fractura, tamano_fractura, ancho_rectangulo, extra_largo_rectangulo);
		extremos_rectangulo		= Extremos_Rectangulo(rectangulo_particulas);
		particulas_fractura		= Potenciales_Particulas(extremos_rectangulo, posicion_particulas, rectangulo_particulas);
		
		//punto_origen_fractura_recortado	= Punto_Origen_Fractura(punto_borde_fractura, tamano_fractura[0] - separacion, angulo_fractura, separacion);
		
		//punto_origen 			= Punto_Origen_Fractura(punto_borde, tamano_fractura[0] - separacion, angulo_fractura, separacion);
		//Util.Imprimir_Vector(extremos_rectangulo);
		//Util.Imprimir_Lista(particulas_fractura);
		//Util.Imprimir_Vector(punto_borde);
		//Util.Imprimir_Vector(punto_origen);
		//System.out.println(posicion_particulas[309][1] + " " + posicion_particulas[310][1]);
		//System.out.println(posicion_particulas[69][1] + " " + posicion_particulas[70][1]);
	}
	
	
	/* 1. crear un metodo donde se cree un rectangulo dado un punto, una longitud, un ancho y un angulo
	 * 
	 * 2. el otro metodo, dependiendo de los maximos y minimos tanto en "x" como en "y", va "filtrando" las
	 *  particulas que esten a una distancia menor a la separacion entre los mismos y las guarda
	 * en un arraylist
	
	 */
		 
	public double[][] Rectangulo(double[] punto_origen, double[] punto_borde, double[] tamano_fractura, double ancho, double sobrante)
	{
		//1. a partir de coordenadas fractura, hallar el origen de la fractura
		
		/*  2.
		* con eso, se tiene los puntos dados por el borde (P_b) de la fractura y el origen (P_o), con esos dos 
		* puntos se "redondea" esta linea de fractura a los nodos mas cercanos que correspondan con esa linea
		* y se guarda en un array usando la posicion_nodos:
		* 
		* 2.1 se crea una especie de rectangulo alrededor de la linea entre los dos puntos, cuyo "exceso"
		* es precisamente la separacion (un poco menor a la mitad de la separacion entre particulas/nodos)
		* 
		* P_b = punto borde,..... P_o = punto origen, cada "-", "_" y "|" tiene una longitud de "separacion"
		* 	"-", y "|"
		* 
		*		 | borde..... 
		*		 |			L_1
		*	   A *------------------------* D
		*		 | *P_b _____________*P_o |
		*	L_4  |						  | L_2
		*	   B * -----------------------* C
		*		 |			L_3
		*		 |
		* 
		* La linea 1 (L_1) son los puntos A y E
		* La linea 2 (L_2) son los puntos D y E
		* La linea 3 (L_3) son los puntos B y D
		* La linea 4 (L_4) son los puntos A y B
		* 
		*/
		
		int    dimension				= Param.dimension;
		double tolerancia 				= 0.01;	//esto es necesario para hacer que el ancho del rectangulo sea un poco mas pequeño que la separacion interparticula para evitar problemas de presicion
		double RAIZ_DOS 				= Param.RAIZ_DOS;
		double PI 						= Param.PI;
		double UNO_ENTRE_RAIZ_DOS 		= (1.0 - tolerancia)/RAIZ_DOS;
		double CUARENTA_CINCO_GRADOS 	= 45.0/180.0*PI;
		double angulo_esquina			= Math.atan(sobrante/ancho);
		double separacion				= Math.sqrt(sobrante*sobrante + ancho*ancho);
		double angulo 					= tamano_fractura[1]/180.0*PI;
		
		double[][] rectangulo = new double[2*dimension][dimension];
		
		/* metodo que a partir de las coordenadas de un punto, una longitud, la separacion interparticula/nodo y
		 * un angulo, crea un rectangulo, la primera componente de esta matriz son los diferentes puntos
		 *  rectangulo[0] -> punto A
		 *  rectangulo[1] -> punto B
		 *  rectangulo[2] -> punto C
		 *  rectangulo[3] -> punto D
		 * la segunda componente indica las componentes "x" y "y"
		 */
	 
		rectangulo[0][0] = punto_borde[0] - separacion*Math.sin(angulo_esquina - angulo);
		rectangulo[0][1] = punto_borde[1] - separacion*Math.cos(angulo_esquina - angulo);
		
		rectangulo[1][0] = punto_borde[0] - separacion*Math.sin(angulo_esquina + angulo);
		rectangulo[1][1] = punto_borde[1] + separacion*Math.cos(angulo_esquina + angulo);

		rectangulo[2][0] = punto_origen[0] - separacion*Math.sin(angulo_esquina + angulo);
		rectangulo[2][1] = punto_origen[1] + separacion*Math.cos(angulo_esquina + angulo);

		rectangulo[3][0] = punto_origen[0] - separacion*Math.sin(angulo_esquina - angulo);
		rectangulo[3][1] = punto_origen[1] - separacion*Math.cos(angulo_esquina - angulo);
		
		return rectangulo;
	}
	
	/* Este metodo retorna los valores maximos y minimos del rectangulo
	*  0 -> valor minimo en x
	*  1 -> valor maximo en x
	*  2 -> valor minimo en y
	*  3 -> valor maximo en y
	*/
	public double[] Extremos_Rectangulo(double[][] rectangulo)
	{
		int    dimension			= Param.dimension;
		
		double[] extremos_rectangulo = new double[2*dimension];

		double a = rectangulo[0][0];
		double b = rectangulo[1][0];
		double c = rectangulo[2][0];
		double d = rectangulo[3][0];
		
		extremos_rectangulo[0] = Math.min(Math.min(a,b),Math.min(c,d));
		extremos_rectangulo[1] = Math.max(Math.max(a,b),Math.max(c,d));
		
		a = rectangulo[0][1];
		b = rectangulo[1][1];
		c = rectangulo[2][1];
		d = rectangulo[3][1];
		
		extremos_rectangulo[2] = Math.min(Math.min(a,b),Math.min(c,d));
		extremos_rectangulo[3] = Math.max(Math.max(a,b),Math.max(c,d));
		
		return extremos_rectangulo;
	}
	
	public List<Integer> Potenciales_Particulas(double[] extremos_rectangulo, double[][] posicion_particulas, double[][] rectangulo)
	{
		
		int longitud	= Param.N;
		int dimension	= Param.dimension;
		
		
		double[] A = Util.Copiar_Vector(rectangulo[0]);
		double[] B = Util.Copiar_Vector(rectangulo[1]);
		double[] D = Util.Copiar_Vector(rectangulo[3]);
		double[] AM = new double[dimension];
		double[] M = new double[dimension];
		double[] AB = Util.Resta_Vectores(A, B);
		double[] AD = Util.Resta_Vectores(A, D);
		double ABAB = Util.Producto_Punto(AB, AB);
		double ADAD = Util.Producto_Punto(AD, AD);

		double AMAB = 0;
		double AMAD = 0;
		

		List<Integer> 	potenciales_particulas = new ArrayList<Integer>();
		
		for(int i = 0; i < longitud; i++)
		{
			M = Util.Copiar_Vector(posicion_particulas[i]);
			AM = Util.Resta_Vectores(A, M);
			AMAB = Util.Producto_Punto(AM, AB);
			AMAD = Util.Producto_Punto(AM, AD);
			if(ABAB > AMAB && AMAB > 0.0 && ADAD > AMAD && AMAD > 0.0)
				potenciales_particulas.add(i);
		}
		
		return potenciales_particulas;
	}
	

	public double[] Punto_Origen_Fractura(double[] coordenadas_fractura, double longitud, double angulo)
	{
		double x, y;
		x = coordenadas_fractura[0] + 0.5*longitud*Math.cos(angulo);
		y = coordenadas_fractura[1] + 0.5*longitud*Math.sin(angulo);
		return new double[] {x, y};
	}

	public double[] Punto_Borde_Fractura(double[] coordenadas_fractura, double longitud, double angulo)
	{
		double x, y;
		x = coordenadas_fractura[0] + 0.5*longitud*Math.cos(angulo + Param.PI);
		y = coordenadas_fractura[1] + 0.5*longitud*Math.sin(angulo + Param.PI);
		return new double[] {x, y};
	}
}

class Propagacion_Fracturas
{
	Parametros 		Param 		= new Parametros();
	Utilidades 		Util 		= new Utilidades();
	
	boolean			hubo_fractura;  
	boolean[][] 	vecinos_hook_nuevo; 
	List<Integer>	particulas_en_borde_nuevo;
	List<Integer>	particulas_fractura_nuevo;
	int[]			borde_de_particula_nuevo;
	double[]		punto_origen_fractura_nuevo;
	double[]		punto_borde_fractura_nuevo;
	double			longitud_fractura_nuevo;
	List<double[]> 	fractura_borde_nuevo;
	List<double[]> 	fractura_origen_nuevo;
	double[]		superficie_nuevo;
	double[][]		vector_normal_nuevo;
	boolean			llego_origen_a_frontera_nuevo;
	boolean			llego_borde_a_frontera_nuevo;

	
	public Propagacion_Fracturas
	(
		int 			t_actual,
		int[][] 		vecinos, 
		double[][] 		estres_prin,
		double[][]		posicion,
		double			separacion,
		double			criterio_de_griffith,
		double[][][]	grad_kernel,
		boolean[][] 	vecinos_hook, 
		List<Integer>	particulas_en_borde,
		List<Integer> 	particulas_fractura, 
		int[]			borde_de_particula,
		double[]		punto_origen_fractura,
		double[]		punto_borde_fractura,
		double			longitud_fractura,
		List<double[]>  fractura_origen,
		List<double[]>  fractura_borde,
		double[]		superficie,
		double[][] 		vector_normal,
		double			longitud_fractura_anterior,
		boolean			llego_origen_a_frontera,
		boolean			llego_borde_a_frontera
	)
	{
						
		//se pasan las referencias a las nuevas variables
		
		hubo_fractura				= false;
		vecinos_hook_nuevo 			= Util.Copiar_Matriz(vecinos_hook);
		particulas_en_borde_nuevo 	= Util.Copiar_Lista(particulas_en_borde);
		particulas_fractura_nuevo	= Util.Copiar_Lista(particulas_fractura);
		borde_de_particula_nuevo	= Util.Copiar_Vector(borde_de_particula);
		punto_origen_fractura_nuevo = Util.Copiar_Vector(punto_origen_fractura);
		punto_borde_fractura_nuevo	= Util.Copiar_Vector(punto_borde_fractura);
		longitud_fractura_nuevo		= longitud_fractura;
		fractura_borde_nuevo		= Util.Copiar_ListaDD(fractura_borde);
		fractura_origen_nuevo		= Util.Copiar_ListaDD(fractura_origen);
		superficie_nuevo			= Util.Copiar_Vector(superficie);
		vector_normal_nuevo			= Util.Copiar_Matriz(vector_normal);


		int 		N 								= Param.N;
		int 		vecina    						= 0;
		int 		cantidad_vecinos 				= 0;
		double[] 	punto_medio 					= new double[Param.dimension];
		double[] 	vector_medio 					= new double[Param.dimension];
		double 		angulo_crecimiento_fractura 	= 0;
		double		angulo_en_que_viene_fractura	= 0;
		double 		longitud_crecimiento_fractura 	= 0;
		double		distancia_media;
		double		modulo;
		double 		tolerancia_distancia			= 1.1;
		double		distancia_maxima				= tolerancia_distancia*Param.RAIZ_DOS*separacion;
		double		distancia_1						= 0;				
		double		distancia_2						= 0;				
		double[]	estres_1 						= Util.Transpuesta(estres_prin)[0];
		double 		tolerancia_estres				= 0.95;
		double		estres_maximo 					= tolerancia_estres*Util.Maximo(estres_1);
		double[] 	Indices_Y_Valores_Maximos_Estres;
		int 		primer_indice_maximo;
		int 		segundo_indice_maximo;
		double 		estres_promedio					= 0;
		double		maxima_diferencia_angulos		= 1./12.*Param.PI;
		boolean		condicion_1						= false;
		boolean		condicion_2						= false;
		
		
		List<Integer> Particulas_Cercanas_Borde 	= new ArrayList<Integer>();
		List<Integer> Particulas_Cercanas_Origen 	= new ArrayList<Integer>();

		for(int i = 0; i < N ; i++)
		{
			distancia_1 = Util.Distancia(posicion[i], punto_origen_fractura);
			distancia_2 = Util.Distancia(posicion[i], punto_borde_fractura);
			
			if(distancia_1 < distancia_maxima && !Particulas_Cercanas_Origen.contains(i))
				Particulas_Cercanas_Origen.add(i);

			if(distancia_2 < distancia_maxima && !Particulas_Cercanas_Borde.contains(i))
				Particulas_Cercanas_Borde.add(i);
		}


		//estres_prin[189][0] = 1;
		//estres_prin[190][0] = 1;

		
		if(!Particulas_Cercanas_Origen.isEmpty() && !llego_origen_a_frontera)
		{
		
			Indices_Y_Valores_Maximos_Estres = Util.Dos_Indices_Y_Elementos_Mayores_En_Array_Con_Lista(Particulas_Cercanas_Origen, Util.Transpuesta(estres_prin)[0]);
			primer_indice_maximo 			 = (int)Indices_Y_Valores_Maximos_Estres[0];
			segundo_indice_maximo 			 = (int)Indices_Y_Valores_Maximos_Estres[1];

			llego_origen_a_frontera_nuevo	 = particulas_en_borde.contains(primer_indice_maximo) && particulas_en_borde.contains(segundo_indice_maximo);

			estres_promedio 				= 0.5*(Indices_Y_Valores_Maximos_Estres[2] + Indices_Y_Valores_Maximos_Estres[3]);
			
			
			condicion_1 = estres_promedio >= criterio_de_griffith && !particulas_fractura.contains(primer_indice_maximo) && !particulas_fractura.contains(segundo_indice_maximo);
			condicion_2 = primer_indice_maximo != -1 && segundo_indice_maximo != -1;
			
			if(condicion_1 && condicion_2)
			{
				
				punto_medio = Util.Promedio_Dos_Vectores(posicion[primer_indice_maximo], posicion[segundo_indice_maximo]);
				
				distancia_1 = Util.Distancia(posicion[primer_indice_maximo], posicion[segundo_indice_maximo]);
				punto_medio = Util.Resta_Vectores(punto_medio, punto_origen_fractura);
				
				angulo_crecimiento_fractura = Util.Angulo_2D(punto_medio);
				angulo_en_que_viene_fractura = Util.Angulo_2D(Util.Resta_Vectores(fractura_origen_nuevo.get(fractura_origen_nuevo.size() - 2), fractura_origen_nuevo.get(fractura_origen_nuevo.size() - 1))); 
				
				if(Math.abs(angulo_crecimiento_fractura - angulo_en_que_viene_fractura) > maxima_diferencia_angulos)
				{
					
					punto_origen_fractura_nuevo[0] = punto_origen_fractura[0] + distancia_1*Math.cos(angulo_crecimiento_fractura);
					punto_origen_fractura_nuevo[1] = punto_origen_fractura[1] + distancia_1*Math.sin(angulo_crecimiento_fractura);
					
					longitud_fractura_nuevo = longitud_fractura + Util.Distancia(punto_origen_fractura_nuevo, punto_origen_fractura);

					superficie_nuevo[primer_indice_maximo]  += separacion;
					superficie_nuevo[segundo_indice_maximo] += separacion;
					
					fractura_origen_nuevo.add(punto_origen_fractura_nuevo);
					
					/*
					if(!particulas_en_borde_nuevo.contains(primer_indice_maximo))
						particulas_en_borde_nuevo.add(primer_indice_maximo);
					if(!particulas_en_borde_nuevo.contains(segundo_indice_maximo))
						particulas_en_borde_nuevo.add(segundo_indice_maximo);
					*/
						
					if(!particulas_fractura_nuevo.contains(primer_indice_maximo))
						particulas_fractura_nuevo.add(primer_indice_maximo);
					if(!particulas_fractura_nuevo.contains(segundo_indice_maximo))
						particulas_fractura_nuevo.add(segundo_indice_maximo);
					
					if(borde_de_particula_nuevo[primer_indice_maximo] == -2)
						borde_de_particula_nuevo[primer_indice_maximo]  = 0;
					
					if(borde_de_particula_nuevo[segundo_indice_maximo] == -2)	
						borde_de_particula_nuevo[segundo_indice_maximo] = 0;
					
					
					
					if(longitud_fractura_anterior == longitud_fractura_nuevo)
						hubo_fractura = false;
					else
						hubo_fractura = true;
						
				}
				else
				{
					hubo_fractura = false;
				}
			}
			else
				hubo_fractura = false; 
		
	
		}
		else
		{
			llego_origen_a_frontera_nuevo = llego_origen_a_frontera;
			hubo_fractura = false;
		}
	
		
		if(!Particulas_Cercanas_Borde.isEmpty() && !llego_borde_a_frontera)
		{
			
			Indices_Y_Valores_Maximos_Estres = Util.Dos_Indices_Y_Elementos_Mayores_En_Array_Con_Lista(Particulas_Cercanas_Borde, Util.Transpuesta(estres_prin)[0]);
			primer_indice_maximo 			 = (int)Indices_Y_Valores_Maximos_Estres[0];
			segundo_indice_maximo 			 = (int)Indices_Y_Valores_Maximos_Estres[1];

			llego_borde_a_frontera_nuevo 	 = particulas_en_borde.contains(primer_indice_maximo) && particulas_en_borde.contains(segundo_indice_maximo);

			estres_promedio 				= 0.5*(Indices_Y_Valores_Maximos_Estres[2] + Indices_Y_Valores_Maximos_Estres[3]);
			
			condicion_1 = estres_promedio >= criterio_de_griffith && !particulas_fractura.contains(primer_indice_maximo) && !particulas_fractura.contains(segundo_indice_maximo);
			condicion_2 = primer_indice_maximo != -1 && segundo_indice_maximo != -1;
			
			if(condicion_1 && condicion_2)
			{
				
				punto_medio = Util.Promedio_Dos_Vectores(posicion[primer_indice_maximo], posicion[segundo_indice_maximo]);
				
				distancia_1 = Util.Distancia(posicion[primer_indice_maximo], posicion[segundo_indice_maximo]);
				punto_medio = Util.Resta_Vectores(punto_medio, punto_borde_fractura);
				
				angulo_crecimiento_fractura = Util.Angulo_2D(punto_medio);
				angulo_en_que_viene_fractura = Util.Angulo_2D(Util.Resta_Vectores(fractura_borde_nuevo.get(fractura_borde_nuevo.size() - 2), fractura_borde_nuevo.get(fractura_borde_nuevo.size() - 1))); 
				
				if(Math.abs(angulo_crecimiento_fractura - angulo_en_que_viene_fractura) > maxima_diferencia_angulos)
				{
					
					punto_borde_fractura_nuevo[0] = punto_borde_fractura[0] + distancia_1*Math.cos(angulo_crecimiento_fractura);
					punto_borde_fractura_nuevo[1] = punto_borde_fractura[1] + distancia_1*Math.sin(angulo_crecimiento_fractura);
					
					longitud_fractura_nuevo = longitud_fractura + Util.Distancia(punto_borde_fractura_nuevo, punto_borde_fractura);

					superficie_nuevo[primer_indice_maximo]  += separacion;
					superficie_nuevo[segundo_indice_maximo] += separacion;
					
					fractura_borde_nuevo.add(punto_borde_fractura_nuevo);
					/*
					if(!particulas_en_borde_nuevo.contains(primer_indice_maximo))
						particulas_en_borde_nuevo.add(primer_indice_maximo);
					if(!particulas_en_borde_nuevo.contains(segundo_indice_maximo))
						particulas_en_borde_nuevo.add(segundo_indice_maximo);
					*/
						
					if(!particulas_fractura_nuevo.contains(primer_indice_maximo))
						particulas_fractura_nuevo.add(primer_indice_maximo);
					if(!particulas_fractura_nuevo.contains(segundo_indice_maximo))
						particulas_fractura_nuevo.add(segundo_indice_maximo);
					
					if(borde_de_particula_nuevo[primer_indice_maximo] == -2)
						borde_de_particula_nuevo[primer_indice_maximo]  = 0;
					
					if(borde_de_particula_nuevo[segundo_indice_maximo] == -2)	
						borde_de_particula_nuevo[segundo_indice_maximo] = 0;
					
					
					
					if(longitud_fractura_anterior == longitud_fractura_nuevo)
						hubo_fractura = false;
					else
						hubo_fractura = true;
						
				}
				else
				{
					hubo_fractura = false;
				}
			}
			else
				hubo_fractura = false;
		
		}
		else
		{
			llego_borde_a_frontera_nuevo = llego_borde_a_frontera;
			hubo_fractura = false;
		}
		
							
		//System.out.println(estres_promedio + " " + Indices_Y_Valores_Maximos_Estres[0] + " " + Indices_Y_Valores_Maximos_Estres[1]);
		
		if(fractura_origen_nuevo.size() != fractura_origen.size() || fractura_borde_nuevo.size() != fractura_borde.size())		
			for(int i = 0; i < N; i++)
			{
				cantidad_vecinos = vecinos[i].length;
				
				for(int j = 0; j < cantidad_vecinos; j++)
				{
					//System.out.println(i  + " " + j + " " + );
					vecina = vecinos[i][j];
					if(vecinos_hook[i][j])
					{
						Interseccion_Lineas Interseccion = new Interseccion_Lineas(posicion[i], posicion[vecina], punto_origen_fractura_nuevo, punto_origen_fractura);
					
						if(Interseccion.se_intersectan)
							vecinos_hook_nuevo[i][j] = false;
						
						Interseccion_Lineas Interseccion2 = new Interseccion_Lineas(posicion[i], posicion[vecina], punto_borde_fractura_nuevo, punto_borde_fractura);
					
						if(Interseccion2.se_intersectan)
							vecinos_hook_nuevo[i][j] = false;
								
						//if(i == 9 &&  vecina == 10) System.out.println(vecinos_hook_nuevo[i][j] );
						
						else
						{
							vector_normal_nuevo[i][0] += -grad_kernel[i][j][0];
							vector_normal_nuevo[i][1] += -grad_kernel[i][j][1];
						}
					}
					else
						vecinos_hook_nuevo[i][j] = false;
				}
				
				modulo = Util.Modulo_Vector(vector_normal_nuevo[i]);
				vector_normal_nuevo[i][0] /= modulo;
				vector_normal_nuevo[i][1] /= modulo;
				
			}
		

		//System.out.println(t_actual + " Si aumento " + fractura.get(fractura.size() - 1)[0] + " " + fractura.get(fractura.size() - 1)[1]);
		//System.out.println(posicion[209][0] + " " + posicion[209][1] + "  " + posicion[210][0] + " " + posicion[210][1]);
	}
	 
	public double[] Direccion_Vector_Normal(double[] posicion_i, double[] posicion_j)
	{
		int dimension = Param.dimension;
		double	modulo = 1;
		double[] vector_centro = new double[dimension];
		
		vector_centro = Util.Producto_Vector_Escalar(0.5, Util.Suma_Vectores(posicion_i, posicion_j));
		vector_centro = Util.Resta_Vectores(posicion_i, vector_centro);
		modulo = Util.Modulo_Vector(vector_centro);
		
		return Util.Producto_Vector_Escalar(1/modulo, vector_centro);
	}
}

class Criterio_Griffith
{
	double criterio_de_griffith;
	
	Parametros Param = new Parametros();
	
	public Criterio_Griffith(double longitud_fractura)
    {
		
		double cociente_tamano_longitud = Param.PI*longitud_fractura/(2.*Param.tamano[0]);
		
		//double cociente_tamano_longitud =  Param.PI*Param.tamano_fractura[0]/(2.*Param.tamano[0]);	
   		//double funcion_auxiliar_griffith	= Math.sqrt(Param.PI*Math.tan(cociente_tamano_longitud)/cociente_tamano_longitud)*(0.752 + 2.02*(2./Param.PI)*cociente_tamano_longitud + 0.37*Math.pow(1 - Math.sin(cociente_tamano_longitud),3.))/Math.cos(cociente_tamano_longitud);	
   		double funcion_auxiliar_griffith	= 1.122  - 0.231*cociente_tamano_longitud + 10.550*Math.pow(cociente_tamano_longitud, 2) - 21.710*Math.pow(cociente_tamano_longitud, 3)  + 30.382*Math.pow(cociente_tamano_longitud, 4);	
		criterio_de_griffith 				= Param.factor_intensidad_estres_critico/(Math.sqrt(longitud_fractura)*funcion_auxiliar_griffith);
	}
}
