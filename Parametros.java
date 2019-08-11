import java.awt.Color;
import java.util.*;

public class Parametros
{


 // 1 Parametros generales
	
	// 1.1. constantes basicas

	static int 	   tmax 			= 0;	// cantidad de iteraciones maximas, si es 0 es solo condiciones iniciales, si es < 0 no hay limite (exepto si tiempo_max_acus = true)
	static int 	   n_sobrante 		= 0;	//cantidad de particulas sobrantes que existe si se coloca UN (1) agujero en el dominio, corregir con este valor a N
	static int 	   nx 				= 120;	//cantidad particulas en "x"
	static int 	   ny 				= 120; //cantidad particulas en "y", para fracturas justo en el centro de la altura en "y", es recomendable usar particulas pares en "y"
	static int 	   N  				= nx*ny - n_sobrante; //cantidad total de particulas
	static boolean elastoestatica	= true; //si es true, aplica condiciones de borde, si es false, no

	
	// 1.2. otras constantes
	
	static double CUATRO_SOBRE_TRES = 1.33333;
	static double PI				= 3.14159;
	static double RAIZ_DOS			= 1.41421;
	static double g 				= 0.0;	//adimensionalizado
    static double angulo_gravedad 	= 90;    //en grados    static double 	delta_t 		= 0.005;        //adimensional
	static int 	  dimension 		= 2;	//por ahora solo implementado para 2D   
	static int	  n_figuras			= 1;	//obsoleto



//2. Parametros del material

	// 2.1 Parametros de escala
	
	static double densidad_real        	  	= 2600; //estos son los parametros reales del material en los cuales se adimensionaliza todo
	static double modulo_compresibi_real  	= 4.62*Math.pow(10,10); //estos son los parametros reales del material en los cuales se adimensionaliza todo
	static double escala_distancia_real		= 1.0/1.0; // en metros
	static double escala_tiempo_real		= Math.sqrt(densidad_real/modulo_compresibi_real)*escala_distancia_real;
	static double escala_energia_real		= modulo_compresibi_real*escala_distancia_real*escala_distancia_real*escala_distancia_real;
	
	//2.2 parametros mecanicos y/o termicos
	static double[] tamano          		= new double[]{ Math.min(nx/(double)ny,1.0),Math.min( ny/(double)nx, 1.0)};	//adimensionalizado, uno (1) por definicion en el eje mas grande
    static double densidad        		  	= 1.0;	//adimensionalizado, uno (1) por definicion
    static double modulo_compresibilidad  	= 1.0;	//adimensionalizado, uno (1) por definicion
	static double modulo_rigidez          	= 2.62*Math.pow(10,10)/modulo_compresibi_real;	//adimensionalizado
	static double vel_sonido              	= Math.sqrt((modulo_compresibilidad + CUATRO_SOBRE_TRES*modulo_rigidez)/densidad); //adimensionalizado
	static double coeficiente_poisson		= 0.5*(3.0*modulo_compresibilidad - 2.0*modulo_rigidez)/(3.0*modulo_compresibilidad + modulo_rigidez);
	static double modulo_young				= 9.0*(modulo_compresibilidad*modulo_rigidez)/(3.0*modulo_compresibilidad + modulo_rigidez);
	static double tension_superficial		= 3.82/escala_energia_real;
	
	    	
    // 2.3 opciones del tiempo
 
	static double 	delta_t			= 0.005;   
	static double   cons_tiempo		= 0.1;	//constante del CFL
	static boolean 	tiempo_variable = true;        //Toma en cuenta del tiempo de fuerza, el CFL y el acustico
    static boolean  tiempo_max_acus = false;		// establece si el tiempo maximo de la simulacion esta dado por el tiempo acustico
    static double   tiempo_acustico = 1.02*Math.max(tamano[0],tamano[1])/vel_sonido;
									 //20;
									 
	//2.4 constantes derivadas de items anteriores
	static double inv_dim 						= 1.0/(double)(dimension + 1);
	static double DOS_RIGIDEZ 					= 2.0*modulo_rigidez;
	static double constante_traza_estres_plano 	= (modulo_compresibilidad - DOS_RIGIDEZ*inv_dim)*DOS_RIGIDEZ/(modulo_compresibilidad + DOS_RIGIDEZ*(1.0 - inv_dim));

	
	
	//2.3 parametros del solido rigido
    static double  	xmax 					= 1.25*Math.max(tamano[0],tamano[1]);    // tamano en x del lienzo, adimensionalizado, se escala a partir del eje mas grande
    static double  	ymax 					= xmax;   // tamano en y del lienzo, adimensionalizado  
    static double[] ubicacion               = new double[]    {xmax/2.0, ymax/2.0};	//adimensionalizado   
	static double	velocidad_angular       = -0.0*vel_sonido;  //adimensionalizado, por ahora no implementado (falta invariancia ante rotaciones)
	static double[] velocidad_polares       = new double[]{0.0*vel_sonido, 0};	//magnitud y angulo de la velocidad inicial de las particulas, adimensionalizado
    static double   angulo_centro_de_masas  = 0;	//el cero es el eje y positivo (que en java apunta hacia abajo)

    
	//2.4 parametros de fractura
	static double   separacion							= 1.0/(double)(nx);
    static double 	cons_hertz 		   					= 36.0*modulo_rigidez*(3*modulo_compresibilidad + modulo_rigidez)/(27*modulo_compresibilidad + 36*modulo_rigidez);
	static double[] tamano_fractura 					= new double[] {10*separacion*tamano[0],0};	//la longitud y el angulo de la linea de fractura
	static double	coordenada_fractura_en_borde		= -0.5*nx*separacion*tamano[0] + 0.5*tamano_fractura[0];
	static double[] coordenadas_fractura 				= new double[] {coordenada_fractura_en_borde, -0*separacion*tamano[1]};	//las coordenadas (x,y) del extremo de la fractura que toca el borde con respecto al centro de masa
	static double[]	ubicacion_orificio					= new double[] {-0.0*tamano[0], 0.00*tamano[1]};	//ESTE VALOR ES CON RESPECTO AL CENTRO DE MASAS
	static double[] semi_ejes_orificio 					= new double[] {0.0*tamano[0], 0.00*tamano[0]};	//tamano del hoyo central en caso de que se desee
    static boolean 	orificio							= new Utilidades().Modulo_Vector(semi_ejes_orificio) != 0;
	static double	angulo_lineas						= 0.0/180.0*PI;
	static int		obstaculo							= -1;	//0 si es una elipse, 1 si son lineas, cualquier otro numero para no tener 
    static double	factor_intensidad_estres_critico	= 10*Math.sqrt((2.0*modulo_young*tension_superficial));

 	double constante_tension 			= Math.pow(2.,11./4.)*Math.sqrt(separacion*modulo_young*tension_superficial*PI)/(PI + 2.);
 	double cociente_tamano_longitud 	= PI*tamano_fractura[0]/(2*tamano[0]);
    //double funcion_auxiliar_griffith	= Math.sqrt(PI*Math.tan(cociente_tamano_longitud)/cociente_tamano_longitud)*(0.752 + 2.02*(2./PI)*cociente_tamano_longitud + 0.37*Math.pow(1 - Math.sin(cociente_tamano_longitud),3.))/Math.cos(cociente_tamano_longitud);	
   	double funcion_auxiliar_griffith	= 1.122  - 0.231*cociente_tamano_longitud + 10.550*Math.pow(cociente_tamano_longitud, 2) - 21.710*Math.pow(cociente_tamano_longitud, 3)  + 30.382*Math.pow(cociente_tamano_longitud, 4);	
	double tension_critica 				= constante_tension/(tamano_fractura[0]*funcion_auxiliar_griffith*funcion_auxiliar_griffith);
	double factor_estres				= 1.5;
	double tension_critica_dimensional 	= factor_estres*tension_critica*modulo_compresibi_real;   

	//2.5 Otros parametros fisicos  
	     
	static double coe_friccion 	  = 0*Math.pow(10,2)*escala_tiempo_real;
    static double energia_elastica_inicial = 0.0;


// 4 Parametros de los modelos computacionales
    	

    
    //4.1 opciones xsph
    
	static double  epsilon_XSPH 	 = 5*Math.pow(10,-2);	//parametro del XSPH
	static int	   xsph_cada_x_pasos = 100;	    //no implementado
	static boolean XSPH 			 = false;	//para las velocidades
	static boolean XSPH_corregido	 = false;
	
	
	// 4.2. Opciones del kernel
	
	static int 	   cons_kern 				= 0;	//0 para Consistencia orden 0, 1 para orden 1, otro numero para sin correccion
	static boolean cons_grad 				= true;	// (REVISAR.) consistencia orden cero para gradiente
	static int 	   metodo_adaptativo 		= 0;	//0 es ninguno, 1 el de monaghan, 2 para el de sigalotti y lopez 2006
	static int 	   kernels 					= 4;	// 0 gausiano (aun no), 1 lucy (aun no), 2 wenland c2 (aun no), 3 spline subico, 4 wendland c4, 5 quintico (aun no), 6 wendland c6 (aun no)
	static double epsilon_kernel_adaptativo = 0.6;
	static double  k_kernel_adaptativo 		= 1.0;
	static double  beta_escala_h 			= 4.0;
	static double  h0 						= Math.pow(N,-1/beta_escala_h);
	static double  parametro_D 				= 2.5;
	static boolean h_separacion				= false;		// Si es "true", fija si la distancia de suavizado se calcula como un porcentaje de la separacion interparticula inicial (cuyo valor numerico se especifica en el archivo Param), 
											//si es "false", el h se calcula como una ley de potencia  de la cantidad de particulas
	static int	   numero_vecinos 			= (int)(nx*ny*PI*h0*h0/tamano[0]/tamano[1]);
	
	//4.3 viscocidad y estres artificial
	
	static double epsilon_estres_artificial 	= 0.3;
	static double n_estres_artificial 			= 4;
    static double epsilon_viscosidad_artificial = 0.1;
    static double alfa_viscocidad_artificial 	= 1.0*Math.sqrt(modulo_compresibi_real)/Math.sqrt(densidad_real)*escala_distancia_real;	
    static double beta_viscocidad_artificial 	= 0.0*modulo_compresibi_real/densidad_real*escala_distancia_real*escala_distancia_real;
    
    
    //4.4  Opciones del modelo de resorte, obsoleto, solo de prueba
    
	static double x_centro 						= xmax/2.0;
	static double y_centro 						= ymax/2.0;
	static double escala_constantes_elasticas	= Math.sqrt(densidad_real)*Math.pow(modulo_compresibi_real,-1.5)/escala_distancia_real;
	static double cons_kx 						= 1.0*Math.pow(10,12)*escala_constantes_elasticas; //adimensionalizado
	static double cons_ky 						= 1.0*Math.pow(10,12)*escala_constantes_elasticas; //adimensionalizado


// 5 Opciones  de las condiciones iniciales y de borde
	
	
	// 5.1 Condiciones iniciales
	
	static int  	tipo_cond_inicial 		= 2;	//obsoleto, 0 para estres constante, 1 para deformacion constante, 2 para elongacion, otro valor para ninguna condicion inicial
	static double[][] deformacion_inicial 	= Deformacion_Inicial(dimension);
	static double[][] estres_inicial 	 	= Estres_Inicial(dimension);
	static double[] cizalla                 = new double[]    {0.0, 0.0};	//paramtros transformacion afin
    static double[] escala                  = new double[]    {1.00, 1.00};	//paramtros transformacion afin
	static int		zona					= 0; //no terminado, establece donde es el "cero" de la deformacion, 0 si es en el centro, 1 para izquierda, 2 para derecha
	
	// 5.2 Condiciones de borde
	/*		
	* Se consideran solo formas rectangulares, cada objeto tiene 4 paredes, cada ciclo de 4 paredes coincide con varias capas de condiciones de 
	* borde sobre esas 4 paredes.
	* 
	* La clasificacion de los bordes es la siguiente:
	*
	* -2 para particulas del bulk - que no estan en el borde-
	*/
		
	static int[] tipo_cond_borde = new int[] 
	{
		1, //indice 0 para particulas que esten en una fractura
		1, //indice  1 para pared izquierda
		1, //indice  2 para pared superior
		1, //indice  3 para pared derecha 
		1 //indice  4 para pared inferior
	};
	/* 
	 *  -1 para particulas a las que no se les haya podido encontrar pared
	* Si el valor de cada elemento es 0, indica que a esa pared no se le aplican condiciones de borde, 1 indica condicion constante de traccion, 
	* 2 indica condicion puntual de traccion, 3 indica condicion constante de dezplazamiento, 4 indica condicion puntual de dezplazamiento
	*/
	static int[] 		orden_condiciones_borde = new int[] {0, 2, 4, 3, 1};
	static int 		   	numero_condiciones_borde = tipo_cond_borde.length;
	/* este vector indica el orden en la cual se ejecutan las condiciones de borde, asi {3, 1, 2, 0, 4} indica que primero la pared 3, luego la 1, etc
	 * por lo que el orden de prioridad mas alto es el ultimo numero (sobreescribe su valor en las particulas de las esquinas)
	 */
	//static double[][] 	vector_traccion 		= Vector_Traccion(numero_condiciones_borde, dimension);
	static double[][]   vector_desplazamiento 	= Vector_Desplazamiento(numero_condiciones_borde, dimension);
	static boolean		condiciones_borde	 	= false;
	
	
	// 5.3 Parametros de las paredes
	
	static int 	  cond_borde			= 1; //para las paredes....0 para inelastico simple, 1 para inelastico mejorado, 2 para pedioricas, 3 para sin delizamiento, otro sin condiciones de borde Solo implementada la opcion 0 por los momentos
	static int    paredes 				= 4;	//implementado solo para paredes cuadradas
	static double pendiente_izquierda 	= 0.0;
	static double pendiente_derecha   	= -0.0;
	static double pendiente_superior  	= 0.0;
	static double pendiente_inferior  	= -0.0;
	static double constante_izquierda 	= 0;
	static double constante_derecha   	= xmax;
	static double constante_superior  	= 0;
	static double constante_inferior  	= ymax ;



// 6 Opciones de visualizacion, aun no clasificado

		//6.1 Opciones del lienzo
    static boolean  graficacion							= true;
    static String 	titulo 								= "SPH";	//titulo del panel de la animacion
	static Color 	color_fondo							= Color.GRAY;	//Color del lienzo 
	static int 		ms_entre_frames 					= 17;	//milisegundos entre cada frame
    static boolean 	escala_color_lateral 				= true;		//escala de colores que aparece a la derecha, por ahora estos colores indica la diferencia relativa entre la densidad
	static boolean 	min_max_totales 					= true;	//opcion de graficacion, si es true, los valores minimos y maximos de la barra de color a la hora de graficar se calcularan de acuerdo al valor mas alto/bajo historico, es decir, a medida de que pase el tiempo se mantendran los valores mas altos o bajos durante el tiempo de simulacion. si es false, estos valores se calcularan como el valor mas alto/bajo de entre las particulas, y se volvera a calcular en cada paso de tiempo
	static int 		parametro_escala_color 				= 0; //no implementado, se espera que sea para cambiar entre grafica de variables en la escala de color, siendo 0 densidad, 1 temperatura, componentes de estres/velocidades, etc
    static Color 	color_marco_inical 					= Color.black;
    static boolean  graficar_marco_inicial 				= false;
    static boolean 	lineas_mallado 						= false;	// REVISAR, hay que cuadrar el mallado con el marco...dibuja la malla de vecinos que tienen ancho kh
	static int 		escala_normalizada					= 600;	//esta escala es solo para visualizacion, ocupa 600 pixeles en la pantalla el lateral mayor, el tamaño del lienzo esta normalizado y no depende del tamaño del dominio computacional
    static int 		escala_lienzo 						= Escala(xmax, ymax, escala_normalizada);
    static double 	ancho_linea_color 					= 0.03;	//indica el ancho de la barra como porcentaje del ancho total de la ventana, fijado a 0.03 -> 3% 
	static double 	ancho_sobre_espacio_lateral_lienzo 	= 1.15;
	static float 	tolerancia_color 					= 0.65f;	//indica el truncamiento del matiz a la hora de graficar la escala de color de la barra y de las particulas, con 0.6 el valor mas bajo es verde y el mas alto es rojo
	static double 	x_max 								= xmax*escala_lienzo;
	static double 	y_max 								= ymax*escala_lienzo;
	static double  	borde_marco 						= 0.1*x_max;
	static int 		x_lienzo 							= (int)(x_max + 2.5*borde_marco + ancho_linea_color*x_max);
	static int 		y_lienzo 							= (int)(y_max + 2.5*borde_marco);

			
	//6.2 Opciones de graficacion de particulas 
	static int 		par_medio_iz 				= (int)(nx*(ny-1)*0.5);
	static int 		par_medio_de 				= (int)(nx*(ny + 1)*0.5) - 1;	//para nx y ny pares
	static int		par_central					= (int)(0.5*(ny)*nx + nx*0.5 );
	static boolean 	graficar_borde 				= true;
    static Color	color_particulas_prueba		= Color.pink;
    static boolean particula_prueba_dominio 	= true;
	static boolean 	modo_particula 				= false;	//indica si se le hace seguimiento a "particula_prueba" (la variable siguiente), mostrando un circulo de tamaño kh y coloreando en negro sus particulas vecinas
	static int 		particula_prueba 	  		//= par_medio_de;
												//= 0*(int)(0.5*(ny)*nx + nx*0.5 );
												= 1952;
	static Color 	color_particulas_borde 		= Color.yellow;
	static boolean  grupo_particulas			= false;
	static boolean  grupo_particulas_dominio	= false;
	//static int[] 	particulas_vistas 			= {1811,1812, 1813, 1814, 1815, 1816, 1817, 1818, 1819, 1820, 1821, 1822, 1823, 1824, 1825, 1826, 1827, 1828, 1829, 1830, 1831, 1832, 1833, 1834, 1835};
	static int[] 	particulas_vistas 			= {1751,1752, 1753, 1754, 1755, 1756, 1757, 1758, 1759, 1760, 1761, 1762, 1763, 1764, 1765, 1766, 1767, 1768, 1769, 1770, 1771, 1772, 1773, 1774, 1775};
	static boolean	graficar_indices_particulas	= false;

	
	//6.3 Opciones de graficacion de los nodos
	static Color	color_nodos							= Color.black;
	static boolean	graficar_nodos						= true;
	static boolean	graficar_indices_nodos				= false;
	static int[]	nodos_en_especifico					//=  {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 0  };
														= {-1}; // -1 para no seleccionar nodo alguno
// 7 Opciones de impresion, aun no clasificado

	static String directorio_base	  	    		= new String(System.getProperty("user.dir"));
	static String directorio_superior				= Directorio_Superior(directorio_base);
	static String directorio_serie_temporal 		= new String(directorio_superior + "serie_temporal/");
	static String nombre_archivo_serie_temp 		= new String("serie_temporal.txt");
	static String serie_temporal 					= new String(directorio_serie_temporal + nombre_archivo_serie_temp);
	static String directorio_serie_espacial 		= new String(directorio_superior + "serie_espacial/");
	static String nombre_archivo_serie_espa 		= new String("serie_espacial");
	static String serie_espacial 					= new String(directorio_serie_espacial + nombre_archivo_serie_espa);
	static int    numero_variables_guardar 			= 43;
	static int    variables_extra_serie_temporal 	= 8;
	static int 	  imprimir_cada_x_pasos		  		= 2;
	static boolean condicion_imprimir				= true;
	static boolean Imprimir_Estres_Maximo_Y_Griffith = true;
	static boolean Imprimir_Info_Fractura			 = true;
	
// 8 Otras opciones, aun no clasificado

	static boolean 	modo_particula_posicion = true;
	static double 	x_par_prueba 			= 0.5*constante_derecha;
	static double 	y_par_prueba 			= 0.5*constante_derecha;
	static boolean 	cerrar_despues_finalizar = false;
	
	static public String Directorio_Superior(String directorio_base)
	{
		String 		separador 			= new String(System.getProperty("file.separator"));
		String[] 	array_directorio 	= directorio_base.split(separador);
		
		String directorio_superior = "";
		for(int i = 0; i < array_directorio.length - 1; i++)
			directorio_superior += array_directorio[i] + separador;
		
		return directorio_superior;
		
	}
    
	static public double[][] Estres_Inicial(int dimension)
    {
		Utilidades Util = new Utilidades();
		
		double[][] estres_inicial = new double[dimension][dimension];
		 

		estres_inicial[0][0] = 0*Math.pow(10,7);

		if(dimension > 1)
		{
			estres_inicial[0][1] = 0;
			estres_inicial[1][0] = estres_inicial[1][0];
			estres_inicial[1][1] = 0;

			if(dimension > 2)
			{
				estres_inicial[0][2] = 0;
				estres_inicial[1][2] = 0;
				estres_inicial[2][0] = estres_inicial[0][2];
				estres_inicial[2][1] = estres_inicial[1][2];
				estres_inicial[2][2] = 0;
			}
		}
		estres_inicial = Util.Producto_Matriz_Escalar(1/modulo_compresibi_real,estres_inicial);
		
	
		return estres_inicial;
	}

	static public double[][] Deformacion_Inicial(int dimension)
    {

		double[][] deformacion_inicial = new double[dimension][dimension];
		

		deformacion_inicial[0][0] = 0*Math.pow(10,-3);

		if(dimension > 1)
		{
			deformacion_inicial[0][1] = 0;
			deformacion_inicial[1][0] = deformacion_inicial[0][1];
			deformacion_inicial[1][1] = 0;

			if(dimension > 2)
			{
				deformacion_inicial[0][2] = 0;
				deformacion_inicial[1][2] = 0;
				deformacion_inicial[2][0] = deformacion_inicial[0][2];
				deformacion_inicial[2][1] = deformacion_inicial[1][2];
				deformacion_inicial[2][2] = 0;
			}
		}
		

		return deformacion_inicial;
	}





	static public double[][] Vector_Desplazamiento(int numero_condiciones_borde, int dimension)
    {
		Utilidades Util = new Utilidades();
		//SE SUPONE QUE EN UNA DIMENSION LA BARRA ESTA ACOSTADA EN EL EJE X
		//la primera componente indica sobre cual pared esta actuando, borde izquierdo = 0, borde superior = 1, borde derecho = 2, borde inferior = 3
		//la segunda indica las componentes espaciales, 0 para "x" y 1 para "y"
		int figuras = 0;	//ahorita hay solo una figura, pero cuando hayan varias y cada una con una condicion de borde diferente, hay que ajustar que cada borde se escale con el modulo de compresibilidad de su figura
		double[][] vector_desplazamiento = new double[numero_condiciones_borde][dimension];
		
		for(int condiciones_borde = 0; condiciones_borde < numero_condiciones_borde; condiciones_borde++)
		{
			if(condiciones_borde == 0)
				vector_desplazamiento[condiciones_borde][0] = 0*Math.pow(10,-2);
			if(condiciones_borde == 2)
				vector_desplazamiento[condiciones_borde][0] = 0*Math.pow(10,-2);

			if(dimension > 1)
			{
				if(condiciones_borde == 0)
					vector_desplazamiento[condiciones_borde][1] = 0*Math.pow(10,6);
				if(condiciones_borde == 2)
					vector_desplazamiento[condiciones_borde][1] = 0*Math.pow(10,6);
				if(condiciones_borde == 1)
				{
					vector_desplazamiento[condiciones_borde][0] = 0*Math.pow(10,-2);
					vector_desplazamiento[condiciones_borde][1] = 0*Math.pow(10,-2);
				}
				if(condiciones_borde == 3)
				{
					vector_desplazamiento[condiciones_borde][0] = 0*Math.pow(10,-2);
					vector_desplazamiento[condiciones_borde][1] = 0*Math.pow(10,-2);
				}

				vector_desplazamiento[condiciones_borde] = Util.Producto_Vector_Escalar(escala_distancia_real,vector_desplazamiento[condiciones_borde]);
			}
		}

		//System.out.println(numero_condiciones_borde); 
		return vector_desplazamiento;	
	}



	static public  final int Escala(double xmax, double ymax, int tamano)
	{   
		if(xmax >= ymax) return (int)(tamano/xmax);
		else return (int)(tamano/ymax);
	
	}

}


