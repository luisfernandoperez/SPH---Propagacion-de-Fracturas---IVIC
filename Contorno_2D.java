public class Contorno_2D
{
    double[] Bordes;
    public Contorno_2D(int figura, double[] tamano, double coordenada)
    {

        if(figura == 0) Bordes = Rectangulo(tamano);
        if(figura == 1) Bordes = Hexagono_Regular(tamano, coordenada);
        if(figura == 2) Bordes = Elipse(tamano, coordenada);


    }


    public double[] Rectangulo(double[] tamano)
    {
        int paredes = 4;
        // 0 para limite izquierdo de x, 1 para limite derecho de x, 2 para limite superior de y, 3 para limite inferior de y
        double[] limites = new double[paredes];
        limites[0] = - 0.5*tamano[0];
        limites[1] = + 0.5*tamano[0];
        limites[2] = - 0.5*tamano[1];
        limites[3] = + 0.5*tamano[1];
        return limites;
    }

    //No terminado
    public double[] Hexagono_Regular(double[] tamano, double coordenada)
    {
        int paredes = 4;
        double UNO_SOBRE_RAIZ_TRES = 0.577350269;
        // 0 para limite izquierdo de x, 1 para limite derecho de x, 2 para limite superior de y, 3 para limite inferior de y
        double[] limites = new double[paredes];

        limites[0] = + UNO_SOBRE_RAIZ_TRES*(Math.abs(coordenada) - tamano[1]);
        limites[1] = + UNO_SOBRE_RAIZ_TRES*(tamano[1] - Math.abs(coordenada));
        limites[2] = - 0.5*tamano[1];
        limites[3] = + 0.5*tamano[1];
        return limites;
    }

    public double[] Elipse(double[] tamano, double coordenada)
    {
        int paredes = 4;
        // 0 para limite izquierdo de x, 1 para limite derecho de x, 2 para limite superior de y, 3 para limite inferior de y
        double[] limites = new double[paredes];
        double termino = tamano[0]/tamano[1]*Math.sqrt(0.25*tamano[1]*tamano[1]- coordenada*coordenada);

        limites[0] = - termino;
        limites[1] = + termino;
        limites[2] = - 0.5*tamano[1];
        limites[3] = + 0.5*tamano[1];
        return limites;
    }

}
