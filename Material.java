//Este metodo calcula las variables intrisecas del material y las asigna a cada particula, posiciones y velocidades CON RESPECTO AL CENTRO DE MASA de cada particula, densidad, velocidad del sonido, separacion y volumen intraparticula, masa de cada particula y energia_elastica . Las constantes del material (compresibilidad y cizalla) son sacados directamente del metodo de Param
// La posicion y velocidad del centro de masa tiene que ser definidas AFUERA de este metodo y se le suman a las posiciones y velocidades de las particulas, es decir, la ubicacion y velocidad de todo el material se define cuando se "invoque"
//Hay que recordar que el material es homogeneo, por lo que todas las particulas tienen los mismos Param intrinsecos
//ESTE METODO NO ESTA ESCRITO EN LA FORMA MAS OPTIMA POSIBLE SINO EN LA MAS SENCILLA DE ENTENDER
//descripcion obsoleta, actualizar
import java.util.Random;
public class Material
{
    double 		volumen;
    double[]	tamano_inicial;
    double[] 	energia_elastica ;
    double   	masas;
    double   	vel_sonido;
    double[]  	separacion;
    double[][] 	posicion_particulas_centro_de_masa;
    double[][] 	velocidad_particulas_centro_de_masa;

    Utilidades 		Util   		= new Utilidades();
    Parametros 		Param 		= new Parametros();
    Otras_Variables Otras_Var 	= new Otras_Variables();


    public Material()
    {

		//Variables basicas
		int N         = Param.N;
        int dimension = Param.dimension;
        
		//Variables macroscopias del material que seran trasladadas a las particulas
        double dens_energia_macro = Param.energia_elastica_inicial;
        double[] tamano           = Param.tamano;  
        
        //Variables principales de salida
        energia_elastica     				 = new double[N];
        separacion 							 = new double[dimension];
        posicion_particulas_centro_de_masa   = new double[N][dimension];
        velocidad_particulas_centro_de_masa  = new double[N][dimension];
        
       	
        //Otras variables
        double[] posicion_centro_de_masa 					= new double[dimension];
        double[][] pos_cent_de_masa_matriz 					= Util.Matriz_Vectores_Iguales(N, posicion_centro_de_masa);
        double[][] posiciones_particulas_centro_geometrico 	= new double[N][dimension];

		//Establecer las variables de las particulas
        separacion 							= Separacion(N, tamano);
        tamano_inicial  					= new double[] {separacion[0], separacion[1], separacion[0]};
        volumen    							= tamano_inicial[0]*tamano_inicial[1]*tamano_inicial[2];
        masas	   							= Param.densidad*volumen;
        energia_elastica     				= Util.Vector_Componentes_Iguales(N, dens_energia_macro);   
        posicion_particulas_centro_de_masa 	= Distribucion_Posiciones(N, separacion, tamano);
        posicion_centro_de_masa 			= Util.Vector_Centro_de_Masas(posicion_particulas_centro_de_masa);
		velocidad_particulas_centro_de_masa = Velocidades_con_Respecto_al_Centro_de_Masas(posicion_particulas_centro_de_masa);	
    }

    public double[] Separacion(int N,  double[] tamano)
    {
        int nx_prima = 0;
        int ny_prima = 0;
        double nx = (double)Param.nx;
        double ny = (double)Param.ny;
        double[] sep = new double[tamano.length + 1];    //0  para sep_eje_x, 2 para sep_eje_y

        sep[0] = tamano[0]/(double)nx;
        sep[1] = sep[0];
        sep[2] = 0.5*(tamano[1] - ny*sep[1]);
        
        return sep;
    }

    public double[][] Distribucion_Posiciones(int N, double[] separacion, double[] tamano)
    {
		boolean 	fuera_dominio 	= false;
        int 		dimension 		= Param.dimension;
        int 		kx 				= 0;
        int 		ky 				= 0;
        int 		k 				= 0;
        int 		kx_ant 			= 0;
        int 		ky_ant 			= 0;
        int 		figura 			= 0;
		int 		n_max 			= Param.nx*Param.ny;
		int 		n_sobrante		= 0;
		double 		py 				= 0;
        double 		px 				= 0;
        double tolerancia 	= 0.01;	//esto es necesario para hacer que el ancho del rectangulo sea un poco mas pequeño que la separacion interparticula para evitar problemas de presicion
		double tol_mas		= 1.0 + tolerancia;
		double tol_menos	= 1.0 - tolerancia;
        double[]	centroide	 	= Param.ubicacion_orificio;
        double[]	semi_ejes 		= Param.semi_ejes_orificio;
        double[]   	pos 			= new double[dimension];
        double[]   	contorno;
        double[][] 	posicion 		= new double[N][dimension];


		double[]   coordenadas_fractura = Param.coordenadas_fractura;
		double[] 	tamano_fractura	  	= Param.tamano_fractura;
		double[] 	punto_origen 		= Punto_Origen_Fractura(coordenadas_fractura, tamano_fractura[0], tamano_fractura[1]);
		
		//rectangulo puro
		double[][] 	rectangulo 			= Rectangulo(punto_origen, coordenadas_fractura, tamano_fractura, 1.0*separacion[0]);
		double[]    extremos			= Extremos_Rectangulo(rectangulo);
		double[] 	V_AB 				= Util.Resta_Vectores(rectangulo[0], rectangulo[1]);
		double[] 	V_BC 				= Util.Resta_Vectores(rectangulo[1], rectangulo[2]);
		double 		T1 					= Util.Producto_Punto(V_AB, V_AB);
		double 		T2 					= Util.Producto_Punto(V_BC, V_BC);
		
		//cuña, rectangulo + triangulo
		double 		pendiente			= Math.tan(Param.angulo_lineas);
		double 		ancho_rectangulo 	= 5.0*separacion[0];
		double 		largo_triangulo  	= 0.5*ancho_rectangulo/pendiente;
		double 		largo_rectangulo 	= tamano_fractura[0] - 1.0*largo_triangulo + 1.5*separacion[0];
		double[] 	punto_origen_cunia  = Punto_Origen_Fractura(coordenadas_fractura, largo_rectangulo, tamano_fractura[1]);
		double[][] 	rectangulo_cunia 	= Rectangulo(punto_origen_cunia, coordenadas_fractura, tamano_fractura, ancho_rectangulo);
		double[]    extremos_cunia		= Extremos_Rectangulo(rectangulo_cunia);
					//extremos_cunia[1]	+= 1*separacion[0];
		double[] 	V_AB_cunia 			= Util.Resta_Vectores(rectangulo_cunia[0], rectangulo_cunia[1]);
		double[] 	V_BC_cunia 			= Util.Resta_Vectores(rectangulo_cunia[1], rectangulo_cunia[2]);
		double 		T1_cunia 			= Util.Producto_Punto(V_AB_cunia, V_AB_cunia);
		double 		T2_cunia 			= Util.Producto_Punto(V_BC_cunia, V_BC_cunia);
		//double[] 	punto_origen_corre  = Punto_Origen_Fractura(coordenadas_fractura, largo_rectangulo + 2.5*separacion[0], tamano_fractura[1]);
		
		for(int i =0; i < n_max; i++)
		{

			py = (ky + 0.5)*separacion[1] - 0.5*tamano[1];
			Contorno_2D  Figura = new Contorno_2D(figura, tamano, py);
			contorno = Figura.Bordes;
			px = contorno[0] + (kx + 0.5)*separacion[0];

			if(px >  contorno[1])
			{
				kx = 0;				
				ky++;
				kx_ant = kx;
				ky_ant = ky;
				pos[0] = contorno[0] + (kx + 0.5)*separacion[0];
				pos[1] = contorno[2] + (ky + 0.5)*separacion[1] + separacion[2];

				kx++;				
			}

			if(px <= contorno[1])
			{
				pos[0] = contorno[0] + (kx + 0.5)*separacion[0];
				pos[1] = contorno[2] + (ky + 0.5)*separacion[1] + separacion[2];
				kx_ant = kx;
				ky_ant = ky;
				kx++;
			}
			
			//System.out.println(i + " " + k);
			
			fuera_dominio = Fuera_Dominio
			(
				k, pendiente,  tol_mas, tol_menos, T1, T2, V_AB, V_BC,
				rectangulo[0], rectangulo[1], T1_cunia, T2_cunia, V_AB_cunia, 
				V_BC_cunia, rectangulo_cunia[0], rectangulo_cunia[1], 
				extremos_cunia, pos,  tamano, centroide, semi_ejes, 
				separacion[0], punto_origen_cunia, tamano_fractura
			);
			
			if(fuera_dominio)
			{
				n_max++;
				i++;
			}
			else
			{
				if(k < N)
				{
					//Random ax = new Random(), ay = new Random();
					ky_ant++;
					//posicion[k][0] = Param.tamano[0]/20*ax.nextGaussian() ;
					//posicion[k][1] = Param.tamano[1]/20*ay.nextGaussian() ;
					posicion[k][0] = contorno[0] + (kx_ant + 0.5)*separacion[0];
					posicion[k][1] = contorno[2] + (ky_ant + 0.5)*separacion[1] + -1*separacion[1];
					k++;
				}
			}
			
					
			
		}

        n_sobrante		= N + n_max - 2*Param.nx*Param.ny;
		//if(Util.Modulo_Vector(semi_ejes) != 0)
		System.out.println(" La cantidad de particulas que sobran del orificio son " + " " + n_sobrante);
		
			for(int i = N - n_sobrante; i < N; i++) 
			{
				posicion[i][0] = -0.65*tamano[0];
				posicion[i][1] = 0.65*tamano[1];
				//System.out.println(i + " " + posicion[i][0]);
				//k++;
			}
			    		
		
			
        return posicion;
        
    }
    
    //este metodo retorna "true" si la particula esta FUERA del dominio definido por los bordes y los obstaculos
    public boolean Fuera_Dominio
    (
		int i,
		double pendiente,
		double tol_mas, 
		double tol_menos, 
		double T1, 
		double T2, 
		double[] V_AB, 
		double[] V_BC, 
		double[] P_A,
		double[] P_B,
		double T1_cunia, 
		double T2_cunia, 
		double[] V_AB_cunia, 
		double[] V_BC_cunia, 
		double[] P_A_cunia,
		double[] P_B_cunia,
		double[] extremos_cunia,
		double[] posicion, 
		double[] tamano, 
		double[] centroide, 
		double[] semi_ejes_orificio, 
		double separacion,
		double[] punto_origen,
		double[] tamano_fractura
	)
    {
		//if(i==51)System.out.println("b");
		
		boolean condicion	= Param.condiciones_borde;
		int obstaculo		= Param.obstaculo;
		
		if(obstaculo == -1)	//sin obstaculos
		{
			condicion = false || posicion[0] > tamano[0] || posicion[1] > tamano[1];
		}
		
		if(obstaculo == 0)	//circulo
		{
			double x			= posicion[0] - centroide[0];
			double y			= posicion[1] - centroide[1];
			double x2			= x*x;
			double y2			= y*y;
			double semi_eje_x2 =  semi_ejes_orificio[0]*semi_ejes_orificio[0];
			double semi_eje_y2 =  semi_ejes_orificio[1]*semi_ejes_orificio[1];

			condicion = x2 < semi_eje_x2*(1.0 - y2/semi_eje_y2) || posicion[0] > tamano[0] || posicion[1] > tamano[1];
		}
		if(obstaculo == 1)	//dos lineas rectas en forma de cuña
		{
			double punto_interseccion_lineas = semi_ejes_orificio[0] - 0.5*tamano[0];
			double recta_1 =   pendiente*(posicion[0] - punto_interseccion_lineas) + centroide[1];
			double recta_2 = - pendiente*(posicion[0] - punto_interseccion_lineas) + centroide[1];
			condicion =  posicion[1] > recta_1 && posicion[1] < recta_2 || posicion[0] > tamano[0] || posicion[1] > tamano[1];

		}
			
		if(obstaculo == 2)	//rectangulo
		{

			double[] V_AD = Util.Resta_Vectores(P_A, posicion);
			double[] V_BD = Util.Resta_Vectores(P_B, posicion);
							

			double T3 = Util.Producto_Punto(V_AB, V_AD);
			double T4 = Util.Producto_Punto(V_BC, V_BD);
						
			condicion = (T1 >= tol_menos*T3 && tol_mas*T3 >= 0.0 && T2 >= tol_menos*T4 && tol_mas*T4 >= 0.0);
		}
		
		if(obstaculo == 3)	//rectangulo con un triangulo en la punta
		{
			/*	
			 * 	 	 __b__________ c
			*	  	|			 |\
			*	  	|			 | \
			*	  a	|			 | /
			*	    |  			 |/
			*	  	 -------------
			* 
			 * a = ancho del rectangulo = ancho del trigangulo
			 * b = largo del rectangulo
			 * c = largo/altura del trigangulo
			 * 
			 * 
			 * */
			boolean condicion_triangulo = false;
			boolean condicion_rectangulo = false;
			
			double recta_1 =   pendiente*(posicion[0] - punto_origen[0]) + punto_origen[1];
			double recta_2 = - pendiente*(posicion[0] - punto_origen[0]) + punto_origen[1];
			condicion_triangulo =  posicion[0] > extremos_cunia[1] &&posicion[1] > recta_1 && posicion[1] < recta_2 || posicion[0] > tamano[0] || posicion[1] > tamano[1];

		
			double[] V_AD_cunia = Util.Resta_Vectores(P_A_cunia, posicion);
			double[] V_BD_cunia = Util.Resta_Vectores(P_B_cunia, posicion);
			double T3 = Util.Producto_Punto(V_AB_cunia, V_AD_cunia);
			double T4 = Util.Producto_Punto(V_BC_cunia, V_BD_cunia);
						
			condicion_rectangulo = (T1_cunia >= tol_menos*T3 && tol_mas*T3 >= 0.0 && T2_cunia >= tol_menos*T4 && tol_mas*T4 >= 0.0);
			
			//if(i==51) System.out.println(punto_origen[0] + " " + recta_1 + " " + recta_2 + " " + extremos_cunia[1] + " " + extremos_cunia[2] + " " + extremos_cunia[3] + " " + condicion_rectangulo);
			
			condicion =   condicion_rectangulo || condicion_triangulo;
		}
		
		return condicion ;
	}

	 
	public double[] Punto_Origen_Fractura(double[] coordenadas_fractura, double longitud, double angulo)
	{
		double x, y;
		x = coordenadas_fractura[0] + longitud*Math.cos(angulo);
		y = coordenadas_fractura[1] + longitud*Math.sin(angulo);
		return new double[] {x, y};
	}
	
	public double[][] Rectangulo(double[] punto_origen, double[] punto_borde, double[] tamano_fractura, double separacion)
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
	 
		rectangulo[0][0] = punto_borde[0] - separacion*UNO_ENTRE_RAIZ_DOS*Math.sin(CUARENTA_CINCO_GRADOS - angulo);
		rectangulo[0][1] = punto_borde[1] - separacion*UNO_ENTRE_RAIZ_DOS*Math.cos(CUARENTA_CINCO_GRADOS - angulo);
		
		rectangulo[1][0] = punto_borde[0] - separacion*UNO_ENTRE_RAIZ_DOS*Math.cos(CUARENTA_CINCO_GRADOS - angulo);
		rectangulo[1][1] = punto_borde[1] + separacion*UNO_ENTRE_RAIZ_DOS*Math.sin(CUARENTA_CINCO_GRADOS - angulo);

		rectangulo[2][0] = punto_origen[0] - separacion*UNO_ENTRE_RAIZ_DOS*Math.cos(CUARENTA_CINCO_GRADOS - angulo);
		rectangulo[2][1] = punto_origen[1] + separacion*UNO_ENTRE_RAIZ_DOS*Math.sin(CUARENTA_CINCO_GRADOS - angulo);

		rectangulo[3][0] = punto_origen[0] - separacion*UNO_ENTRE_RAIZ_DOS*Math.sin(CUARENTA_CINCO_GRADOS - angulo);
		rectangulo[3][1] = punto_origen[1] - separacion*UNO_ENTRE_RAIZ_DOS*Math.cos(CUARENTA_CINCO_GRADOS - angulo);
		
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
	

    public double[][] Velocidades_con_Respecto_al_Centro_de_Masas(double[][] posiciones)
    {
		int 		N 				  = posiciones.length;
        int 		dimension  		  = posiciones[0].length;
        double 		angulo 			  = 0;
        double 		velocidad_angular = 0;
        double 		modulo_distancia  = 0;
        double[] 	posicion_rotada	= new double[dimension];
        double[][] 	velocidad 		  = new double[N][dimension];
        
        for(int i =0; i < N; i++)
        {
			angulo 				= Util.Angulo_2D(posiciones[i]);
			velocidad_angular 	= Param.velocidad_angular;
			modulo_distancia	= Util.Modulo_Vector(posiciones[i]);
			posicion_rotada = Util.Producto_Matriz_Vector(Util.Derivada_Matriz_Rotacion_2D(angulo), posiciones[i]);
			velocidad[i] = Util.Producto_Vector_Escalar(velocidad_angular, posicion_rotada);


        }
        return velocidad;
    }

}
