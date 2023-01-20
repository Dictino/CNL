# CNL
Herramientas del curso de [CNL](http://portal.uned.es/portal/page?_pageid=93,70656202&_dad=portal&_schema=PORTAL&idAsignatura=31104178&idTitulacion=310401) del [Máster en Ingeniería de Sistemas y Control UNED/UCM](https://cv4.ucm.es/moodle/course/view.php?id=4056)

## Instalación
Para ejecutar los ejemplos lo primero que necesitamos es instalar Julia, para ello lo mejor es ir a https://julialang.org/ y visitar su página de [descargas](https://julialang.org/downloads/) y descargar a última versión estable.

La instalación es trivial, en windows o mac basta con ejecutar el fichero descargado.

En Linux basta con descomprimir el fichero tar en una carpeta cualquiera de vuestra elección y después ejecutar julia que se encuentra en la subcarpeta /bin/julia. (Por comodidad es conveniente crear un enlace simbólico desde la carpeta donde se encuentra julia a una carpeta que esté en el path ej. ```sudo ln -s $(pwd)/julia /bin/julia```)

A continuación, desgargad los ficheros en vuestro ordenador y hecho esto abrid julia, ir a la carpeta donde están descargados los ficheros:


```julia
julia> cd("C:\\ruta\\en\\mi\\ordenador")

```

NOTA: podéis usar la barra espaciadora para ir autocompletando la dirección o copiar y pegar desde el explorador y luego cambiar los / por //. Otra opción es crear un enlace directo a julia en esa carpeta o en una direccion del path y abrir dicho enlace desde el navegador.


Esto **tardará un buen rato** la primera vez ya que tiene que decargar e instalar todos los paquetes que usará el curso (Pluto,  Plots, OrdinaryDiffEq ModelingToolkit, SymEngine...) así que aprovechad para ver algún tema de CNL ;)

Cuando termine se abrirá pluto en el navegador:

![pluto](https://github.com/Dictino/CNL/blob/main/Im%C3%A1genes/pluto.png?raw=true)

Desde ahí se puede seleccionar el notebook que se desea abrir pinchando sobre la barra, buscando el fichero y pulsando Open.

Si queréis probarlo **sin instalar nada** en vuestro ordenador (no es la opción recomendada, pero es muy útil si queréis verlo en una tableta, teléfono o en un ordenador en el que no queréis instalar software) usad el siguiente enlace

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/Dictino/CNL/HEAD?urlpath=pluto)

Ojo que binder **no es persistente**, si hacéis un cambio y queréis guardarlo no uséis el botón de guardar sino el de exportar:

![pluto](https://github.com/Dictino/CNL/blob/main/Im%C3%A1genes/Boton_exportar.png?raw=true)

¡A disfrutarlo!

Una guía para conocer la sintaxis básica de Julia es el siguiente [tutorial en castellano](https://hedero.webs.upv.es/julia-basico/), o también [estos en inglés](https://julialang.org/learning/tutorials/) (el primero "From zero to Julia!" merece bastante la pena).

Un recurso fundamental es la [documentación](https://docs.julialang.org/en/v1/) que es muy completa y exausiva y si se quiere concoer el lenguaje el manual es de obligada lectura, para los que dominéis otros lenguajes puede ser interesante la siguiente [sección](https://docs.julialang.org/en/v1/manual/noteworthy-differences/).
 
