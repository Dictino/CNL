md"""

# CNL
Herramientas del curso de [CNL](http://portal.uned.es/portal/page?_pageid=93,70656202&_dad=portal&_schema=PORTAL&idAsignatura=31104178&idTitulacion=310401) del [Máster en Ingeniería de Sistemas y Control UNED/UCM](https://cv4.ucm.es/moodle/course/view.php?id=4056)

## Instalación
Para ejecutar los ejemplos lo primero que necesitamos es instalar Julia, para ello lo mejor es ir a https://julialang.org/downloads/ y visitar su página de [descargas](https://julialang.org/downloads/) y descacargar a última versión estable. La instalación es trivial, basta con ejecutar el fichero descargado o en el caso de Linux basta con descomprimir el fichero tar en una carpeta cualquiera de vuestra elección y después ejecutar julia que se encuenta en la subcarpeta /bin/julia. (Para comodidad es convenietne crear un simbólico desde la carpeta donde se encuentra julia a una carpeta que esté en el path ej. ```sudo ln -s $(pwd)/julia /bin/julia```)

A continuacion desde el terminal de julia:

```julia
julia>] add Pluto PlutoUI Plots OrdinaryDiffEq ModelingToolkit SymEngine LaTeXStrings

```

Esto tardará un rato (sobre todo con mala conexión a internet y si se usa windows) así que aprovechad para ver algún tema de CNL.

Cuando acabe podeis lanzar pluto:

```julia
julia> using Pluto
julia> Pluto.run()
```

Y se abrirá pluto en el navegador, desde ahí se puede seleccionar el notebook que se desea abrir (poniendo por ejemplo la url de github de dicho notebook o descargándolo en el ordenador y usando dihcha ruta).

![pluto](https://github.com/Dictino/CNL/blob/main/Im%C3%A1genes/pluto.png?raw=true)

¡A disfrutarlo!

"""
