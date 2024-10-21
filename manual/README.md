# CORA manual

This folder contains the source files for the CORA manual.

## How to build the manual

Please see the respective MAKEFILE:

    - MAKEFILE_LINUX
    - MAKEFILE_WINDOWS

**Important:** The CORA manual contains many tikz images.
Thus, building the manual may take a while the first time you run it (<30min). 
It is important to enable shell escape (`--shell-escape`) 
to not run into memory issues. During the first execution, all tikz figures
are saved as pdf, making subsequent executions much faster (<20s).

Additionally, you might want to increase the TeX main memory size.
On windows, execute in the cmd:

    initexmf --edit-config-file=pdflatex

and add the following line to the opened file:

    main_memory=12000000

The value `12000000` is arbitrary, but it works for me.

Finally, you may need to execute the following command in the cmd:

	initexmf --dump=pdflatex

### Building chapter-wise

Each chapter in the manual is included in the main document (`CoraManual.tex`)
via the `\include` command. This allows an easy compilation of a 
single chapter, which shortens compilation time.
To build only a single chapter, uncomment the `\includeonly` command near
the top of the main document and specify the respective chapter.

<hr style="height: 1px;">

<img src="../app/images/coraLogo_readme.svg"/>