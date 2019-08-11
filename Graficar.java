import java.awt.*;
import java.awt.geom.Ellipse2D;
import java.awt.geom.Rectangle2D;
import java.awt.geom.Line2D;
import java.util.*;
import java.util.List;

public class Graficar
{
	Utilidades Util  = new Utilidades();
	Parametros Param = new Parametros();
	
		//variables del dominio computacional y las particulas
	double xmax = Param.xmax;
	double ymax = Param.ymax;
	double r;


	//variables del lienzo
	int escala 			= Param.escala_lienzo;
	double  borde_marco = Param.borde_marco;
	int x_lienzo 		= Param.x_lienzo;
	int y_lienzo 		= Param.y_lienzo;
	int dimension		= Param.dimension;

	double[] lineas_lienzo_i = new double[dimension];
	double[] lineas_lienzo_f = new double[dimension];
    
	//variables de la graficacion de los observables
	boolean 	min_max_totales 		 = Param.min_max_totales;
	int[] 		particula_medio 		 = new int[1];
	int[] 		particula_medio_superior = new int[1];
	int[] 		particula_medio_inferior = new int[1];
	double[] 	observable;
	float 		maximo;
	float 		minimo;
	float 		min_temp;
	float 		max_temp;
	int 		tamano_letra;
	boolean		graficar_indices_particulas = Param.graficar_indices_particulas;
	boolean 	graficar_indices_nodos		= Param.graficar_indices_nodos;
			
	//variables de graficaciondel panel lateral
	int 	posicion_letras_x 					= (int)(0.85*x_lienzo);
	int 	posicion_letras_y 					= (int)(0.95*y_lienzo);
	float 	matiz_linea 						= 1.0f;
	double 	ancho_linea_color 					= Param.ancho_linea_color*Param.x_max;
	boolean escala_color_lateral 				= Param.escala_color_lateral;
	double 	ancho_sobre_espacio_lateral_lienzo 	= Param.ancho_sobre_espacio_lateral_lienzo;

	//variables de la graficacion del mallado de vecinos

	double 	 posiciones_ajustada_de_lineas_borde = 0;	//0 para mostrar solo celdas reales, 2 para mostrar celdas reales y fantasmas
	int 	 factor_celdas 						 = 1;
	Color    color_nodos						 = Param.color_nodos;
	double	 radio_nodos						 = 0;
	boolean  graficar_nodos						 = Param.graficar_nodos;
	int[]	 nodos_en_especifico				 = Param.nodos_en_especifico;

	//variables de la graficacion de las particulas
	float 	tolerancia_color 						= Param.tolerancia_color;
	float 	matiz 									= 1f;
	int 	r_escala 								= (int)(escala*r);
    int 	particula_prueba 						= Param.particula_prueba;
    boolean lineas_mallado 							= Param.lineas_mallado;
    boolean modo_particula 							= Param.modo_particula;
    boolean grupo_particulas						= Param.grupo_particulas;
    boolean grupo_particulas_dominio				= Param.grupo_particulas_dominio;
	Color 	color 									= Color.black;
	Color 	color_matiz 							= Color.black;
	int[] 	particulas_vistas 						= Param.particulas_vistas;
	double 	soporte_kernel;
	int 	n_graficacion 							= Param.N;
	boolean graficar_borde 							= Param.graficar_borde;
	boolean graficar_marco_inicial					= Param.graficar_marco_inicial;
	Color 	color_particulas_borde 					= Param.color_particulas_borde;
	Color	color_particulas_prueba					= Param.color_particulas_prueba;
	boolean 	particula_prueba_dominio				= Param.particula_prueba_dominio;
	Color 	color_marco_inical 						= Param.color_marco_inical;
	int    	n_figuras 								= Param.n_figuras;
	double 	radio_particulas 						= 2.*escala*r;
	double 	factor_correccion_graficar_particulas 	= borde_marco - radio_particulas*0.5;
	double 	factor_correccion_graficar_rectangulo_x = borde_marco;
	double 	factor_correccion_graficar_rectangulo_y = borde_marco;
	double 	x_centro_rectangulo;
	double 	y_centro_rectangulo;
	double 	tamano_x_rectangulo;
	double 	tamano_y_rectangulo;	
	double pos_x_particulas;
	double pos_y_particulas;
	double xo, xf, yo, yf;
	
	//Otras variables
	boolean condicion_borde = Condicion_Borde();
	double 	k_kernel 		= new Kernel_y_Derivadas().K_Kernel();

	
	Graphics2D   g_2D;
	int 		 t_actual;
	int 		 tmax;
	double 		 tiempo_real;
	double[][] 	 posicion;
	double[][] 	 elongacion;	  	  
	double[] 	 posicion_centro_de_masa;  	  
	double[][] 	 posicion_particulas_centro_de_masa;
	double[][] 	 velocidad;
	double[] 	 velocidad_centro_de_masa;	  
	double[][] 	 velocidad_particulas_centro_de_masa;	  
	double   	 angulo_centro_de_masas;	  
	double[][] 	 aceleracion;  
	double[] 	 energia_elastica ;
	double[] 	 tasa_energia_elastica ;		  
	double[][][] tensor_deformacion;	  
	double[][][] tensor_estres;  
	double[]  	 angulo_estres;		  
	double[][] 	 estres_prin;	  		
	double[]  	 h;
	double 		 ancho_de_celda;
	int[][]		 vecinos;
	List<Integer>	particulas_en_borde;
	boolean[][]  vecinos_hook;
	double[][] 	 posicion_nodos;
	int[] 		 nodos_borde;
	int[]		 nodos_entre_particula;
	double[] 	 densidad;
	int[]		borde_de_particula;
	List<double[]> fractura_origen;
	List<double[]> fractura_borde;
		
	public Graficar
	(
		Graphics2D   g_2D,
		int 		 t_actual,
		int 		 tmax,
		double 		 tiempo_real,
		double[][] 	 posicion,
		double[][] 	 elongacion,	  	  
		double[] 	 posicion_centro_de_masa,  	  
		double[][] 	 posicion_particulas_centro_de_masa,
		double[][] 	 velocidad,
		double[] 	 velocidad_centro_de_masa,	  
		double[][] 	 velocidad_particulas_centro_de_masa,	  
		double[][] 	 aceleracion,  
		double[] 	 energia_elastica ,
		double[] 	 tasa_energia_elastica ,		  
		double[][][] tensor_deformacion,	  
		double[][][] tensor_estres,  
		double[][] 	 estres_prin,	  		
		double[]  	 h,
		double 		 ancho_de_celda,
		int[][]		 vecinos,
		List<Integer>	particulas_en_borde,
		double 		r,
		boolean[][] vecinos_hook,
		double[][] 	posicion_nodos,
		double[] 	densidad,
		int[] 		borde_de_particula,
		List<double[]> fractura_origen,
		List<double[]> fractura_borde
	)
	{
		this.g_2D								= g_2D;
		this.t_actual 							= t_actual;
		this.tmax 								= tmax;
		this.tiempo_real 						= tiempo_real;
		this.posicion 							= posicion;
		this.elongacion 						= elongacion;	  
		this.posicion_centro_de_masa 			= posicion_centro_de_masa;  	  
		this.posicion_particulas_centro_de_masa = posicion_particulas_centro_de_masa;
		this.velocidad 							= velocidad;
		this.velocidad_centro_de_masa 			= velocidad_centro_de_masa;	  
		this.velocidad_particulas_centro_de_masa = velocidad_particulas_centro_de_masa;  
		this.aceleracion 						= aceleracion; 
		this.energia_elastica 		 					= energia_elastica ;
		this.tasa_energia_elastica 		 				= tasa_energia_elastica ;		   
		this.tensor_deformacion 				= tensor_deformacion;  
		this.tensor_estres 						= tensor_estres;
		this.estres_prin 						= estres_prin; 		
		this.h 									= h;
		this.ancho_de_celda 					= ancho_de_celda;
		this.vecinos 							= vecinos;
		this.particulas_en_borde 				= particulas_en_borde;
		this.r 									= r;
		this.vecinos_hook						= vecinos_hook;
		this.posicion_nodos						= posicion_nodos;
		this.densidad							= densidad;
		this.borde_de_particula					= borde_de_particula;
		this.fractura_origen					= fractura_origen;
		this.fractura_borde						= fractura_borde;
		
		
		if(t_actual < tmax + 1 | tmax < 0)
		{
				
			
			//System.out.println(tensor_estres.length);
			observable = Util.Copiar_Vector
			(
				//grad_kernel[0]
				//Util.Suma_Vectores(Util.Transpuesta(estres_prin)[0], Util.Transpuesta(estres_prin)[1])
				Util.Transpuesta(estres_prin)[0]
				//Util.Inv_indices_Matriz(tensor_estres)[0][0]
				//Util.Inv_indices_Matriz(tensor_estres)[1][1]
				//Util.Inv_indices_Matriz(tensor_deformacion)[0][0]
				//Util.Inv_indices_Matriz(tensor_estres)[0][1]
				//Util.Inv_indices_Matriz(tensor_deformacion)[0][0]
				//Util.Suma_Vectores(Util.Inv_indices_Matriz(tensor_estres)[0][0], Util.Inv_indices_Matriz(tensor_estres)[1][1])
				//Util.Resta_Vectores(Util.Inv_indices_Matriz(tensor_estres)[0][0],Util.Inv_indices_Matriz(tensor_estres)[1][1])
				//Util.Vector_de_Traza_R2(tensor_estres)
                // Util.Transpuesta(velocidad)[1]
				//Util.Transpuesta(posicion)[0]
                // Util.Transpuesta(aceleracion)[1]
                 //Util.Transpuesta(aceleracion)[0]
                //Util.Transpuesta(elongacion)[0]
               //Util.Transpuesta(elongacion)[1]
                //Util.Error_Array(funcion, funcion_SPH, funcion)
                //funcion
                //funcion_SPH
                //suma_kernel
				//vel_sonido
				//h
				//Util.Inv_indices_Matriz(tensor_deformacion)[1][1]
				//masas
				//suma_kernel
				//Util.Transpuesta(vecinos)[2]
				//densidad
			);
			


			
			//System.out.println(observable[0] + " " + tensor_deformacion[0][0][0]+ " " + t_actual);
            
            if( min_max_totales)
            {
                if(t_actual ==0)
                {
                    maximo = (float)Util.Maximo(Util.Copiar_Vector_Hasta_N(observable, n_graficacion));
                    minimo = (float)Util.Minimo(Util.Copiar_Vector_Hasta_N(observable, n_graficacion));
                }
				  //System.out.println(minimo + " " + Util.Minimo(Util.Copiar_Vector_Hasta_N(observable, n_graficacion)));
				if(t_actual > 0)
                {
                maximo = (float)Util.Maximo_2_numeros(Util.Maximo(Util.Copiar_Vector_Hasta_N(observable, n_graficacion)), maximo);
                minimo = (float)Util.Minimo_2_numeros(Util.Minimo(Util.Copiar_Vector_Hasta_N(observable, n_graficacion)), minimo);
				}

            }

            if(!min_max_totales)
            {
                maximo = (float)Util.Maximo(Util.Copiar_Vector_Hasta_N(observable, n_graficacion));
                minimo = (float)Util.Minimo(Util.Copiar_Vector_Hasta_N(observable, n_graficacion));
            }
         }
	}
	
	public void Lienzo()
	{
		g_2D.setColor(Param.color_fondo);	//color del fondo
		g_2D.fillRect(0, 0, x_lienzo, y_lienzo);
		// Dibuja el dominio computacional
		if(graficar_marco_inicial)
		{
			lineas_lienzo_i[0] = borde_marco/escala ;
			lineas_lienzo_i[1] = borde_marco/escala ;
			lineas_lienzo_f[0] = xmax +  borde_marco/escala;
			lineas_lienzo_f[1] = ymax +  borde_marco/escala;
			g_2D.setColor(color);
			g_2D.draw(new Line2D.Double((Paredes(lineas_lienzo_i)[0])*escala + borde_marco,lineas_lienzo_i[1]*escala, (Paredes(lineas_lienzo_f)[0])*escala + borde_marco, lineas_lienzo_f[1]*escala));
			g_2D.draw(new Line2D.Double((Paredes(lineas_lienzo_i)[1])*escala + borde_marco,lineas_lienzo_i[1]*escala, (Paredes(lineas_lienzo_f)[1])*escala + borde_marco, lineas_lienzo_f[1]*escala));
			g_2D.draw(new Line2D.Double(lineas_lienzo_i[0]*escala,(Paredes(lineas_lienzo_i)[2])*escala + borde_marco, lineas_lienzo_f[0]*escala, (Paredes(lineas_lienzo_f)[2])*escala + borde_marco));
			g_2D.draw(new Line2D.Double(lineas_lienzo_i[0]*escala,(Paredes(lineas_lienzo_i)[3])*escala + borde_marco, lineas_lienzo_f[0]*escala, (Paredes(lineas_lienzo_f)[3])*escala + borde_marco));
			g_2D.draw(new Line2D.Double(borde_marco,borde_marco + ymax*escala/2,borde_marco + xmax*escala, borde_marco + ymax*escala/2));
			g_2D.draw(new Line2D.Double(borde_marco + xmax*escala/2,borde_marco,borde_marco + xmax*escala/2, borde_marco + ymax*escala));
		}
		
		if(graficar_marco_inicial)
		{	
			g_2D.setColor(color_marco_inical);
			

			tamano_x_rectangulo = escala*Param.tamano[0];
			tamano_y_rectangulo = escala*Param.tamano[1];	
			x_centro_rectangulo = escala*posicion_centro_de_masa[0] + factor_correccion_graficar_rectangulo_x - tamano_x_rectangulo*0.5;
			y_centro_rectangulo = escala*posicion_centro_de_masa[1] + factor_correccion_graficar_rectangulo_y - tamano_y_rectangulo*0.5;
			g_2D.draw(new Rectangle2D.Double(x_centro_rectangulo, y_centro_rectangulo, tamano_x_rectangulo, tamano_y_rectangulo ));
		}
	}
	
	public void Graficar_Particulas()
	{
		//grafica las particulas
		
		if(t_actual < tmax + 1 | tmax < 0)
		{
     			
			radio_particulas = escala*r*1;
			factor_correccion_graficar_particulas = borde_marco - radio_particulas*0.5;

			for(int i =0; i<observable.length; i++)
			{
			
				if(maximo != minimo) matiz =  (float)(tolerancia_color/(maximo-minimo)*(observable[i]-minimo)) + (1-tolerancia_color);
				else matiz = 1f-tolerancia_color/2f;
				double prom = 0.5*(maximo + minimo);
				color_matiz = Color.getHSBColor(matiz, 1, 1);
				/*
				if(observable[i] < prom)
					color_matiz = Color.getHSBColor(0f, 1, 1);
				if(observable[i] > prom)
					color_matiz = Color.getHSBColor(0.5f, 1, 1);
				*/
				xo = escala*posicion[i][0] + factor_correccion_graficar_particulas;
				yo = escala*posicion[i][1] + factor_correccion_graficar_particulas;

				g_2D.setColor(color_matiz);
				g_2D.fill(new Ellipse2D.Double(xo, yo, radio_particulas, radio_particulas) );			
			 	
			
				if(graficar_borde)
				if(borde_de_particula[i] > - 2)
				{	
					//System.out.println(i);
					g_2D.setColor(Brillo(color_matiz, 0.5));
					xo = escala*posicion[i][0] + factor_correccion_graficar_particulas;
					yo = escala*posicion[i][1] + factor_correccion_graficar_particulas;
					g_2D.fill(new Ellipse2D.Double(xo, yo, radio_particulas, radio_particulas) );				
				}
				
			}
			
			
		}

	}
	
	public void Graficar_Numeros_Particulas()
	{
		if(t_actual < tmax + 1 | tmax < 0 && graficar_indices_particulas)
		{
     			
			radio_particulas = escala*r*1;
			factor_correccion_graficar_particulas = borde_marco - radio_particulas*0.5;
			
			for(int i =0; i<posicion.length; i++)
			{
				tamano_letra = 11;
				xo = escala*posicion[i][0] + factor_correccion_graficar_particulas + 0.4*radio_particulas;
				yo = escala*posicion[i][1] + factor_correccion_graficar_particulas + 0.75*radio_particulas;
				if(i >= 100)
				{
					tamano_letra = 10;
					xo = escala*posicion[i][0] + factor_correccion_graficar_particulas + 0.25*radio_particulas;
				}
				if(i >= 10000)
				{
					tamano_letra = 9;
					xo = escala*posicion[i][0] + factor_correccion_graficar_particulas + 0.25*radio_particulas;
				}

				g_2D.setColor(Color.white);
				g_2D.setFont(new Font("TimesRoman", Font.PLAIN, tamano_letra)); 
				g_2D.drawString(String.valueOf(i),(int)(xo), (int)(yo));
			}
		}
	}
	
	public void Graficar_Nodos()
	{
		//grafica los nodos o puntos entre las particulas
		
		if(t_actual < tmax + 1 | tmax < 0 && graficar_nodos)
		{
			//Util.Imprimir_Matriz(posicion_nodos);
			
			radio_nodos = escala*r*0.3;
			soporte_kernel = 3*radio_nodos;
			factor_correccion_graficar_particulas = borde_marco - radio_nodos*0.5;
			
			//System.out.println(posicion_nodos.length + " a" );
			
			for(int i =0; i< posicion_nodos.length; i++)
			{
				g_2D.setColor(color_nodos);
				//if(Util.Indice_De_Elemento(i, nodos_borde) != - 1) /*|| Util.Indice_De_Elemento(i, nodos_entre_particula) != - 1 ) */
				//	g_2D.setColor(color_particulas_borde);
				
				xo = escala*posicion_nodos[i][0] + factor_correccion_graficar_particulas;
				yo = escala*posicion_nodos[i][1] + factor_correccion_graficar_particulas;
				
				//System.out.println(i + " " + xo + " " + yo);
				
				g_2D.fill(new Ellipse2D.Double(xo, yo, radio_nodos, radio_nodos) );
			}
							
			for(int j: nodos_en_especifico)	
			{	
				if(j > -1)
				{
					xo = escala*posicion_nodos[j][0] + borde_marco -soporte_kernel;
					yo = escala*posicion_nodos[j][1] + borde_marco  -soporte_kernel;
					g_2D.draw(new Ellipse2D.Double(xo, yo, 2*soporte_kernel, 2*soporte_kernel) );
				}
			}
				
			
		}

	}
	
	public void Graficar_Fractura()
	{
		g_2D.setColor(Color.black);
		for(int i = 0; i < fractura_borde.size() - 1; i++)
		{
			xo = escala*fractura_borde.get(i)[0] + borde_marco;
			yo = escala*fractura_borde.get(i)[1] + borde_marco;
			xf = escala*fractura_borde.get(i + 1)[0] + borde_marco;
			yf = escala*fractura_borde.get(i + 1)[1] + borde_marco;
			g_2D.draw(new Line2D.Double(xo,yo,xf, yf));
		}
		for(int i = 0; i < fractura_origen.size() - 1; i++)
		{
			xo = escala*fractura_origen.get(i)[0] + borde_marco;
			yo = escala*fractura_origen.get(i)[1] + borde_marco;
			xf = escala*fractura_origen.get(i + 1)[0] + borde_marco;
			yf = escala*fractura_origen.get(i + 1)[1] + borde_marco;
			g_2D.draw(new Line2D.Double(xo,yo,xf, yf));
		}
		
		xo = escala*fractura_origen.get(0)[0] + borde_marco;
		yo = escala*fractura_origen.get(0)[1] + borde_marco;
		xf = escala*fractura_borde.get(0)[0] + borde_marco;
		yf = escala*fractura_borde.get(0)[1] + borde_marco;
		g_2D.draw(new Line2D.Double(xo,yo,xf, yf));
	}
	
	public void Graficar_Numeros_Nodos()
	{
		if(t_actual < tmax + 1 | tmax < 0 && graficar_indices_nodos)
		{
     			
			radio_particulas = escala*r*0.3;
			factor_correccion_graficar_particulas = borde_marco - radio_particulas*0.5;
			
			for(int i =0; i<posicion_nodos.length; i++)
			{
				tamano_letra = 11;
				xo = escala*posicion_nodos[i][0] + factor_correccion_graficar_particulas;
				yo = escala*posicion_nodos[i][1] + factor_correccion_graficar_particulas + 2*radio_particulas;
				if(i >= 100)
				{
					tamano_letra = 10;
					xo = escala*posicion_nodos[i][0] + factor_correccion_graficar_particulas;
				}
				if(i >= 10000)
				{
					tamano_letra = 9;
					xo = escala*posicion_nodos[i][0] + factor_correccion_graficar_particulas + radio_particulas;
				}

				g_2D.setColor(Color.red);
				g_2D.setFont(new Font("TimesRoman", Font.PLAIN, tamano_letra)); 
				g_2D.drawString(String.valueOf(i),(int)(xo), (int)(yo));
			}
		}
	}
	
	public void Auto_Ajuste()
	{
		// auto ajuste del marco del lienzo
		if(ancho_de_celda > 0.3*x_lienzo)
		{
			borde_marco = factor_celdas*(0.3*x_lienzo)/2.0;
			posiciones_ajustada_de_lineas_borde = factor_celdas*ancho_de_celda - factor_celdas*0.3*x_lienzo/2.0 ;
		}
	}

	public void Barra_Lateral()
	{
		if(escala_color_lateral)
		{
			for(int m = 0; m < y_lienzo; m++)
			{
				matiz_linea = ((y_lienzo-tolerancia_color*m)/(y_lienzo));
			    color_matiz = Color.getHSBColor(matiz_linea, 1, 1);
			    g_2D.setColor(color_matiz);
			    g_2D.draw(new Line2D.Double(x_lienzo - ancho_linea_color,m, x_lienzo,m));
			}
			g_2D.setColor(color);
			g_2D.drawString(Float.toString(maximo), posicion_letras_x, y_lienzo -posicion_letras_y);
			g_2D.drawString(Float.toString(minimo), posicion_letras_x, posicion_letras_y);
		}
	}

	public void Mallado()
	{
		// lineas del mallado
		if(lineas_mallado)
		{
			g_2D.setColor(color);
			for(int lx = 0; lx < x_lienzo; lx++){g_2D.draw(new Line2D.Double(ancho_de_celda*lx-posiciones_ajustada_de_lineas_borde,0,ancho_de_celda*lx-posiciones_ajustada_de_lineas_borde,y_lienzo));}
			for(int ly = 0; ly < y_lienzo; ly++){g_2D.draw(new Line2D.Double(0,ancho_de_celda*ly-posiciones_ajustada_de_lineas_borde,x_lienzo,ancho_de_celda*ly-posiciones_ajustada_de_lineas_borde));}

		}
	}



	public void Graficar_Particulas_En_Especifico()
	{
		int vecino;
		if(modo_particula)
		{
			radio_particulas = escala*r;
			factor_correccion_graficar_particulas = borde_marco - radio_particulas*0.5;
			
			int k = Param.particula_prueba;

			xo = escala*posicion[k][0] + factor_correccion_graficar_particulas;
			yo = escala*posicion[k][1] + factor_correccion_graficar_particulas;
			g_2D.setColor(Brillo(color_particulas_prueba, 1.2));
			g_2D.fill(new Ellipse2D.Double(xo, yo, radio_particulas, radio_particulas) );
				
			for(int j = 0; j < vecinos[k].length; j++)
			{				
				vecino = vecinos[k][j];	
				
				if( vecino != k && particula_prueba_dominio && vecinos_hook[k][j])
				//if( vecino != k && particula_prueba_dominio)
				{
					g_2D.setColor(color);
					xo = escala*posicion[vecino][0] + factor_correccion_graficar_particulas;
					yo = escala*posicion[vecino][1] + factor_correccion_graficar_particulas;
					g_2D.fill(new Ellipse2D.Double(xo, yo, radio_particulas, radio_particulas) );
				}
			}
			soporte_kernel = escala*k_kernel*h[k];
			g_2D.setColor(Color.white);
			xo = escala*posicion[k][0] + borde_marco -soporte_kernel;
			yo = escala*posicion[k][1] + borde_marco  -soporte_kernel;
			g_2D.draw(new Ellipse2D.Double(xo, yo, 2.0*soporte_kernel, 2.0*soporte_kernel) );
			
		}
		
		if(grupo_particulas)
		{
		
			for(int k:Param.particulas_vistas)
			{	
				//System.out.println(k + " " + vecinos[k].length);
			    for(int j = 0; j < vecinos[k].length; j++)
			    {
					vecino = vecinos[k][j];

					//if( k == vecino && vecinos_hook[k][j])
					if( k == vecino)
					{
						
						g_2D.setColor(color_particulas_borde);
						
						if(Util.Pertenece_A_La_Lista(vecino, vecinos[Param.particula_prueba]))
							 g_2D.setColor(Brillo(color_particulas_borde, 0.5));
						if((k==Param.particula_prueba))
							 g_2D.setColor(Brillo(color_particulas_prueba, 0.5));
							
					//	System.out.println(k + " ");	
						xo = escala*posicion[vecino][0] + factor_correccion_graficar_particulas;
						yo = escala*posicion[vecino][1] + factor_correccion_graficar_particulas;
						g_2D.fill(new Ellipse2D.Double(xo, yo, radio_particulas, radio_particulas) );
					} 

					//if(vecino !=k && grupo_particulas_dominio && vecinos_hook[k][j])
					if(vecino !=k && grupo_particulas_dominio && vecinos_hook[k][j])
					{
					    g_2D.setColor(color);
					    xo = escala*posicion[vecino][0] + factor_correccion_graficar_particulas;
					    yo = escala*posicion[vecino][1] + factor_correccion_graficar_particulas;
						g_2D.fill(new Ellipse2D.Double(xo, yo, radio_particulas, radio_particulas) );

					}
					
			    }
			    
			    if(grupo_particulas_dominio)
				{
					soporte_kernel = escala*k_kernel*h[k];
					g_2D.setColor(color);
					xo = escala*posicion[k][0] + borde_marco -soporte_kernel;
					yo = escala*posicion[k][1] + borde_marco  -soporte_kernel;
					g_2D.draw(new Ellipse2D.Double(xo, yo, 2.0*soporte_kernel, 2.0*soporte_kernel) );
				}
			}
		}
	}



	
	public void Graficar_Puntos()
	{
		//grafica los nodos o puntos entre las particulas
		
		if(t_actual < tmax + 1 | tmax < 0)
		{
			int o, o_1;
			int longitud = nodos_entre_particula.length  - 1;
			
			g_2D.setColor(color_particulas_borde);
			
			for(int i = 0;  i < longitud;i++)
			{
				o = nodos_entre_particula[i];
				o_1 = nodos_entre_particula[i + 1];
				
				xo = escala*posicion_nodos[o][0] + borde_marco;
				yo = escala*posicion_nodos[o][1] + borde_marco;
				xf = escala*posicion_nodos[o_1][0] + borde_marco;
				yf = escala*posicion_nodos[o_1][1] + borde_marco;
				
				g_2D.draw(new Line2D.Double(xo, yo, xf, yf));

			}
				
			
		}

	}
	
	
	private static Color Brillo(Color c, double scale) 
	{
		int r = Math.min(255, (int) (c.getRed() * scale));
		int g = Math.min(255, (int) (c.getGreen() * scale));
		int b = Math.min(255, (int) (c.getBlue() * scale));
		return new Color(r,g,b);
	}

	public int[] Particula_Medio(double[][] posicion_particulas_centro_de_masa, double[][] separacion)
	{
		int n_figuras = posicion_particulas_centro_de_masa.length;
		int N = posicion_particulas_centro_de_masa.length;
		int[]  particula_medio = new int[n_figuras];
		
		for(int figura = 0; figura < n_figuras; figura++)
		{
			for(int i = 0; i < N; i++)
			{
				if( Math.abs(posicion_particulas_centro_de_masa[i][0]) < separacion[0][0] && Math.abs(posicion_particulas_centro_de_masa[i][1]) < separacion[0][1]) 
					particula_medio[figura] = i;
			}
		} 
		return particula_medio;
	}
	
		public int[] Particula_Medio_Superior(double[][] posicion_particulas_centro_de_masa, double[][] separacion, double[][] tamano)
	{
		int n_figuras = posicion_particulas_centro_de_masa.length;
		int N = posicion_particulas_centro_de_masa.length;
		int[]  particula_medio_superior = new int[n_figuras];
		
		for(int figura = 0; figura < n_figuras; figura++)
		{
			for(int i = 0; i < N; i++)
			{
				if( Math.abs(posicion_particulas_centro_de_masa[i][0]) < separacion[0][0] && Math.abs(posicion_particulas_centro_de_masa[i][1]) < separacion[0][1]) 
					particula_medio_superior[figura] = i;
			}
		} 
		return particula_medio_superior;
	}
	
		public int[] Particula_Medio_Inferior(double[][] posicion_particulas_centro_de_masa, double[][] separacion)
	{
		int n_figuras = posicion_particulas_centro_de_masa.length;
		int N = posicion_particulas_centro_de_masa.length;
		int[]  particula_medio_inferior = new int[n_figuras];
		
		for(int figura = 0; figura < n_figuras; figura++)
		{
			for(int i = 0; i < N; i++)
			{
				if( Math.abs(posicion_particulas_centro_de_masa[i][0]) < separacion[0][0] && Math.abs(posicion_particulas_centro_de_masa[i][1]) < separacion[0][1]) 
					particula_medio_inferior[figura] = i;
			}
		} 
		return particula_medio_inferior;
	}
	
	

	
		// Paredes computacionales, no la del sistema fisico
    public final double[] Paredes(double[] posicion)
    {
        double[] paredes = new double[Param.paredes];
        // 0: pared izquierda, 1: pared derecha, 2: pared superior, 3: pared inferior,
        paredes[0] = Param.pendiente_izquierda*posicion[1] + Param.constante_izquierda;
        paredes[1] = Param.pendiente_derecha*posicion[1]   + Param.constante_derecha;
        paredes[2] = Param.pendiente_superior*posicion[0]  + Param.constante_superior;
        paredes[3] = Param.pendiente_inferior*posicion[0]  + Param.constante_inferior;
        return paredes;
    }

	public boolean Condicion_Borde()
	{
		boolean condicion_borde = false;
		int[] 	tipo_cond_borde = Param.tipo_cond_borde;
		for(int i = 0; i < tipo_cond_borde.length; i++)
		{
			condicion_borde = tipo_cond_borde[i] != 0;
			if(condicion_borde)
				i = tipo_cond_borde.length; 
		}
		return condicion_borde;
	}
}
