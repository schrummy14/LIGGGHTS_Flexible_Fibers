import os
import glob 

def buildCommands():
    lines = """


<!DOCTYPE html>
<!--[if IE 8]><html class="no-js lt-ie9" lang="en" > <![endif]-->
<!--[if gt IE 8]><!--> <html class="no-js" lang="en" > <!--<![endif]-->
<head>
<meta charset="utf-8">

<meta name="viewport" content="width=device-width, initial-scale=1.0">

<title>4. Commands &mdash; LIGGGHTS v3.X documentation</title>











    





    <link rel="stylesheet" href="_static/css/theme.css" type="text/css" />





    <link rel="top" title="LIGGGHTS v3.X documentation" href="index.html"/>
        <link rel="next" title="5. Contact models" href="Section_gran_models.html"/>
        <link rel="prev" title="3. Input Script" href="Section_input_script.html"/> 


<script src="_static/js/modernizr.min.js"></script>

</head>

<body class="wy-body-for-nav" role="document">


<div class="wy-grid-for-nav">

    
    <nav data-toggle="wy-nav-shift" class="wy-nav-side">
    <div class="wy-side-scroll">
        <div class="wy-side-nav-search">
        

        
            <a href="Manual.html" class="icon icon-home"> LIGGGHTS
        

        
        </a>

        
            
            
            <div class="version">
                v3.X
            </div>
            
        

        
<div role="search">
<form id="rtd-search-form" class="wy-form" action="search.html" method="get">
    <input type="text" name="q" placeholder="Search docs" />
    <input type="hidden" name="check_keywords" value="yes" />
    <input type="hidden" name="area" value="default" />
</form>
</div>

        
        </div>

        <div class="wy-menu wy-menu-vertical" data-spy="affix" role="navigation" aria-label="main navigation">
        
            
            
            
            
            
            <ul class="current">
<li class="toctree-l1"><a class="reference internal" href="Section_intro.html">1. Introduction</a></li>
<li class="toctree-l1"><a class="reference internal" href="Section_start.html">2. Getting Started</a></li>
<li class="toctree-l1"><a class="reference internal" href="Section_input_script.html">3. Input Script</a></li>
<li class="toctree-l1 current"><a class="current reference internal" href="">4. Commands</a><ul>
<li class="toctree-l2"><a class="reference internal" href="#list-of-all-commands">4.1. List of all commands</a></li>
<li class="toctree-l2"><a class="reference internal" href="#bond-style-potentials">4.2. bond_style potentials</a></li>
<li class="toctree-l2"><a class="reference internal" href="#compute-styles">4.3. compute styles</a></li>
<li class="toctree-l2"><a class="reference internal" href="#dump-styles">4.4. dump styles</a></li>
<li class="toctree-l2"><a class="reference internal" href="#fix-styles">4.5. fix styles</a></li>
<li class="toctree-l2"><a class="reference internal" href="#pair-style-potentials">4.6. pair_style potentials</a></li>
</ul>
</li>
<li class="toctree-l1"><a class="reference internal" href="Section_gran_models.html">5. Contact models</a></li>
<li class="toctree-l1"><a class="reference internal" href="Section_mesh_modules.html">6. Mesh modules</a></li>
<li class="toctree-l1"><a class="reference internal" href="Section_packages.html">7. Packages</a></li>
<li class="toctree-l1"><a class="reference internal" href="Section_howto.html">8. How-to discussions</a></li>
<li class="toctree-l1"><a class="reference internal" href="Section_modify.html">9. Modifying &amp; extending LIGGGHTS(R)-PUBLIC</a></li>
<li class="toctree-l1"><a class="reference internal" href="Section_python.html">10. Python interface to LIGGGHTS(R)-PUBLIC</a></li>
<li class="toctree-l1"><a class="reference internal" href="Section_errors.html">11. Errors</a></li>
</ul>

            
        
        </div>
    </div>
    </nav>

    <section data-toggle="wy-nav-shift" class="wy-nav-content-wrap">

    
    <nav class="wy-nav-top" role="navigation" aria-label="top navigation">
        
        <i data-toggle="wy-nav-top" class="fa fa-bars"></i>
        <a href="Manual.html">LIGGGHTS</a>
        
    </nav>


    
    <div class="wy-nav-content">
        <div class="rst-content">
        















<div role="navigation" aria-label="breadcrumbs navigation">

<ul class="wy-breadcrumbs">
    
    <li><a href="Manual.html">Docs</a> &raquo;</li>
        
    <li>4. Commands</li>
    
    
    <li class="wy-breadcrumbs-aside">
        
            
            <a href="_sources/Section_commands.txt" rel="nofollow"> View page source</a>
        
        <a href="http://www.cfdem.com"> Website</a>
        
            <a href="Section_commands.html#comm" rel="nofollow"> Commands</a>
            
        
        
    </li>
    
</ul>


<hr/>

    <div class="rst-footer-buttons" role="navigation" aria-label="footer navigation">
    
        <a href="Section_gran_models.html" class="btn btn-neutral float-right" title="5. Contact models" accesskey="n">Next <span class="fa fa-arrow-circle-right"></span></a>
    
    
        <a href="Section_input_script.html" class="btn btn-neutral" title="3. Input Script" accesskey="p"><span class="fa fa-arrow-circle-left"></span> Previous</a>
    
    </div>

</div>
        <div role="main" class="document" itemscope="itemscope" itemtype="http://schema.org/Article">
        <div itemprop="articleBody">
            
<div class="section" id="commands">
<h1>4. Commands<a class="headerlink" href="#commands" title="Permalink to this headline">¶</a></h1>
<p>This section describes how a LIGGGHTS(R)-PUBLIC input script is formatted and the
input script commands used to define a LIGGGHTS(R)-PUBLIC simulation.</p>
<div class="contents local topic" id="contents">
<ul class="simple">
<li><a class="reference internal" href="#list-of-all-commands" id="id1">List of all commands</a></li>
<li><a class="reference internal" href="#bond-style-potentials" id="id2">bond_style potentials</a></li>
<li><a class="reference internal" href="#compute-styles" id="id3">compute styles</a></li>
<li><a class="reference internal" href="#dump-styles" id="id4">dump styles</a></li>
<li><a class="reference internal" href="#fix-styles" id="id5">fix styles</a></li>
<li><a class="reference internal" href="#pair-style-potentials" id="id6">pair_style potentials</a></li>
</ul>
</div>
<div class="section" id="list-of-all-commands">
<span id="comm"></span><span id="cmd-1"></span><h2><a class="toc-backref" href="#id1">4.1. List of all commands</a><a class="headerlink" href="#list-of-all-commands" title="Permalink to this headline">¶</a></h2>
<p>This section lists all LIGGGHTS commands alphabetically, with a separate
listing below of styles within certain commands. Note
that some style options for some commands are part of
packages, which means they cannot be used unless the package was
included when LAMMPS was built.  Not all packages are included in a
default build.  These dependencies are listed as Restrictions
in the command&#8217;s documentation.</p>
<table border="1" class="docutils">
<colgroup>
<col width="16%" />
<col width="16%" />
<col width="16%" />
<col width="18%" />
<col width="18%" />
<col width="18%" />
</colgroup>
<tbody valign="top">
<tr class="row-odd"><td><a class="reference internal" href="atom_modify.html"><em>atom_modify</em></a></td>
<td><a class="reference internal" href="atom_style.html"><em>atom_style</em></a></td>
<td><a class="reference internal" href="bond_coeff.html"><em>bond_coeff</em></a></td>
<td><a class="reference internal" href="bond_style.html"><em>bond_style</em></a></td>
<td><a class="reference internal" href="boundary.html"><em>boundary</em></a></td>
<td><a class="reference internal" href="box.html"><em>box</em></a></td>
</tr>
<tr class="row-even"><td><a class="reference internal" href="change_box.html"><em>change_box</em></a></td>
<td><a class="reference internal" href="clear.html"><em>clear</em></a></td>
<td><a class="reference internal" href="communicate.html"><em>communicate</em></a></td>
<td><a class="reference internal" href="compute.html"><em>compute</em></a></td>
<td><a class="reference internal" href="compute_modify.html"><em>compute_modify</em></a></td>
<td><a class="reference internal" href="create_atoms.html"><em>create_atoms</em></a></td>
</tr>
<tr class="row-odd"><td><a class="reference internal" href="create_box.html"><em>create_box</em></a></td>
<td><a class="reference internal" href="delete_atoms.html"><em>delete_atoms</em></a></td>
<td><a class="reference internal" href="delete_bonds.html"><em>delete_bonds</em></a></td>
<td><a class="reference internal" href="dielectric.html"><em>dielectric</em></a></td>
<td><a class="reference internal" href="dimension.html"><em>dimension</em></a></td>
<td><a class="reference internal" href="displace_atoms.html"><em>displace_atoms</em></a></td>
</tr>
<tr class="row-even"><td><a class="reference internal" href="dump.html"><em>dump</em></a></td>
<td><a class="reference internal" href="dump_modify.html"><em>dump_modify</em></a></td>
<td><a class="reference internal" href="echo.html"><em>echo</em></a></td>
<td><a class="reference internal" href="fix.html"><em>fix</em></a></td>
<td><a class="reference internal" href="fix_modify.html"><em>fix_modify</em></a></td>
<td><a class="reference internal" href="group.html"><em>group</em></a></td>
</tr>
<tr class="row-odd"><td><a class="reference internal" href="if.html"><em>if</em></a></td>
<td><a class="reference internal" href="include.html"><em>include</em></a></td>
<td><a class="reference internal" href="info.html"><em>info</em></a></td>
<td><a class="reference internal" href="jump.html"><em>jump</em></a></td>
<td><a class="reference internal" href="label.html"><em>label</em></a></td>
<td><a class="reference internal" href="lattice.html"><em>lattice</em></a></td>
</tr>
<tr class="row-even"><td><a class="reference internal" href="log.html"><em>log</em></a></td>
<td><a class="reference internal" href="mass.html"><em>mass</em></a></td>
<td><a class="reference internal" href="neigh_modify.html"><em>neigh_modify</em></a></td>
<td><a class="reference internal" href="neigh_modify.html"><em>neigh_settings</em></a></td>
<td><a class="reference internal" href="neighbor.html"><em>neighbor</em></a></td>
<td><a class="reference internal" href="neighbor.html"><em>neighbor_skin</em></a></td>
</tr>
<tr class="row-odd"><td><a class="reference internal" href="newton.html"><em>newton</em></a></td>
<td><a class="reference internal" href="next.html"><em>next</em></a></td>
<td><a class="reference internal" href="orient.html"><em>orient</em></a></td>
<td><a class="reference internal" href="origin.html"><em>origin</em></a></td>
<td><a class="reference internal" href="pair_coeff.html"><em>pair_coeff</em></a></td>
<td><a class="reference internal" href="pair_style.html"><em>pair_style</em></a></td>
</tr>
<tr class="row-even"><td><a class="reference internal" href="partition.html"><em>partition</em></a></td>
<td><a class="reference internal" href="print.html"><em>print</em></a></td>
<td><a class="reference internal" href="processors.html"><em>processors</em></a></td>
<td><a class="reference internal" href="quit.html"><em>quit</em></a></td>
<td><a class="reference internal" href="read_data.html"><em>read_data</em></a></td>
<td><a class="reference internal" href="read_dump.html"><em>read_dump</em></a></td>
</tr>
<tr class="row-odd"><td><a class="reference internal" href="read_restart.html"><em>read_restart</em></a></td>
<td><a class="reference internal" href="region.html"><em>region</em></a></td>
<td><a class="reference internal" href="replicate.html"><em>replicate</em></a></td>
<td><a class="reference internal" href="reset_timestep.html"><em>reset_timestep</em></a></td>
<td><a class="reference internal" href="restart.html"><em>restart</em></a></td>
<td><a class="reference internal" href="run.html"><em>run</em></a></td>
</tr>
<tr class="row-even"><td><a class="reference internal" href="run_style.html"><em>run_style</em></a></td>
<td><a class="reference internal" href="set.html"><em>set</em></a></td>
<td><a class="reference internal" href="shell.html"><em>shell</em></a></td>
<td><a class="reference internal" href="thermo.html"><em>thermo</em></a></td>
<td><a class="reference internal" href="thermo_modify.html"><em>thermo_modify</em></a></td>
<td><a class="reference internal" href="thermo_style.html"><em>thermo_style</em></a></td>
</tr>
<tr class="row-odd"><td><a class="reference internal" href="timestep.html"><em>timestep</em></a></td>
<td><a class="reference internal" href="uncompute.html"><em>uncompute</em></a></td>
<td><a class="reference internal" href="undump.html"><em>undump</em></a></td>
<td><a class="reference internal" href="unfix.html"><em>unfix</em></a></td>
<td><a class="reference internal" href="units.html"><em>units</em></a></td>
<td><a class="reference internal" href="variable.html"><em>variable</em></a></td>
</tr>
<tr class="row-even"><td><a class="reference internal" href="velocity.html"><em>velocity</em></a></td>
<td><a class="reference internal" href="write_data.html"><em>write_data</em></a></td>
<td><a class="reference internal" href="write_dump.html"><em>write_dump</em></a></td>
<td><a class="reference internal" href="write_restart.html"><em>write_restart</em></a></td>
<td>&nbsp;</td>
<td>&nbsp;</td>
</tr>
</tbody>
</table>
</div>"""
    
    lines = lines.split('\n')
    return lines

def buildBonds():
    lines = """<hr class="docutils" />
<div class="section" id="bond-style-potentials">
<h2><a class="toc-backref" href="#id2">4.2. bond_style potentials</a><a class="headerlink" href="#bond-style-potentials" title="Permalink to this headline">¶</a></h2>
<p>See the <a class="reference internal" href="bond_style.html"><em>bond_style</em></a> command for an overview of bond
potentials.  Click on the style itself for a full description:</p>
<table border="1" class="docutils">
<colgroup>"""
    lines = lines.split('\n')
    bondTypes = sorted(glob.glob('bond_*.html'))
    bondTypes.remove("bond_coeff.html")
    bondTypes.remove("bond_style.html")
    bondTypes.remove("bond_none.html")
    bondTypes.remove("bond_hybrid.html")
    numBondTypes = len(bondTypes)+2
    colWidth = int(100/numBondTypes)
    lines.append('<col width="%i%%" >' % (100 - (numBondTypes-1)*colWidth))

    for _ in range(1,numBondTypes):
        lines.append('<col width="%i%%" >' % colWidth)
    
    lines = lines + ["""</colgroup>""", """<tbody valign="top">""", """<tr class="row-odd">"""]
    for bondType in bondTypes:
        lines.append("""<td><a class="reference internal" href="%s"><em>%s</em></a></td>""" % (bondType, bondType[5:-5]))
    lines.append("""<td><a class="reference internal" href="bond_hybrid.html"><em>hybrid</em></a></td>""")
    lines.append("""<td><a class="reference internal" href="bond_none.html"><em>none</em></a></td>""")

    lines.append("</tr>")
    lines.append("</tbody>")
    lines.append("</table>")
    lines.append("</div>")

    return lines


def buildComputes():
    lines = """<div class="section" id="compute-styles">
<h2><a class="toc-backref" href="#id3">4.3. compute styles</a><a class="headerlink" href="#compute-styles" title="Permalink to this headline">¶</a></h2>
<p>See the <a class="reference internal" href="compute.html"><em>compute</em></a> command for one-line descriptions of
each style or click on the style itself for a full description:</p>
<table border="1" class="docutils">
<colgroup>
<col width="34%" />
<col width="33%" />
<col width="33%" />
</colgroup>
<tbody valign="top">"""
    lines = lines.split('\n')

    computeTypes = sorted(glob.glob('compute_*.html'))
    k = 0
    kk = 0
    for computeType in computeTypes:
        if k % 3 == 0:
            if k > 0:
                lines.append("</tr>")
            if kk % 3 == 0 or kk % 3 == 2:
                lines.append("""<tr class="row-odd">""")
            else:
                lines.append("""<tr class="row-even">""")
            kk += 1
        k += 1
        computeName = computeType[8:].replace("_","/").replace(".html","")
        curStr = """<td><a class="reference internal" href="%s"><em>%s</em></a></td>""" % (computeType, computeName)
        lines.append(curStr)
    kRemain = 3 - (k%3)
    if kRemain < 3:
        for k in range(kRemain):
            lines.append("""<td>&nbsp;</td>""")
    lines.append("""</tr>""")
    lines.append("""</tbody>""")
    lines.append("""</table>""")
    lines.append("""</div>""")

    return lines

def buildDumps():
    lines = """<div class="section" id="dump-styles">
<h2><a class="toc-backref" href="#id4">4.4. dump styles</a><a class="headerlink" href="#dump-styles" title="Permalink to this headline">¶</a></h2>
<p>See the <a class="reference internal" href="dump.html"><em>dump</em></a> command for one-line descriptions
of each style or click on the style itself for a full description:</p>
<table border="1" class="docutils">
<colgroup>
<col width="27%" />
<col width="20%" />
<col width="33%" />
<col width="20%" />
</colgroup>
<tbody valign="top">
<tr class="row-odd"><td><a class="reference internal" href="dump_custom_vtk.html"><em>custom/vtk</em></a></td>
<td><a class="reference internal" href="dump_image.html"><em>image</em></a></td>
<td><a class="reference internal" href="dump_local_gran_vtk.html"><em>local/gran/vtk</em></a></td>
<td><a class="reference internal" href="dump_image.html"><em>movie</em></a></td>
</tr>
</tbody>
</table>
</div>"""
    lines = lines.split('\n')
    return lines

def buildFixes():
    lines = """<div class="section" id="fix-styles">
<h2><a class="toc-backref" href="#id5">4.5. fix styles</a><a class="headerlink" href="#fix-styles" title="Permalink to this headline">¶</a></h2>
<p>See the <a class="reference internal" href="fix.html"><em>fix</em></a> command for one-line descriptions
of each style or click on the style itself for a full description:</p>
<table border="1" class="docutils">
<colgroup>
<col width="34%" />
<col width="33%" />
<col width="33%" />
</colgroup>
<tbody valign="top">"""
    lines = lines.split('\n')

    fixTypes = sorted(glob.glob('fix_*.html'))
    fixTypes.remove("fix_modify.html")

    k = 0
    kk = 0
    for fixType in fixTypes:
        if k % 3 == 0:
            if k > 0:
                lines.append("</tr>")
            if kk % 3 == 0 or k % 3 == 2:
                lines.append("""<tr class="row-odd">""")
            else:
                lines.append("""<tr class="row-even">""")
            kk += 1
        k += 1
        fixName = fixType[4:].replace("_","/").replace(".html","")
        curStr = """<td><a class="reference internal" href="%s"><em>%s</em></a></td>""" % (fixType, fixName)
        lines.append(curStr)
    kRemain = 3 - (k%3)
    if kRemain > 0:
        for k in range(kRemain):
            lines.append("""<td>&nbsp;</td>""")
    lines.append("""</tr>""")
    lines.append("""</tbody>""")
    lines.append("""</table>""")
    lines.append("""</div>""")

    return lines

def buildPairStyles():
    lines = """<div class="section" id="pair-style-potentials">
<h2><a class="toc-backref" href="#id6">4.6. pair_style potentials</a><a class="headerlink" href="#pair-style-potentials" title="Permalink to this headline">¶</a></h2>
<p>See the <a class="reference internal" href="pair_style.html"><em>pair_style</em></a> command for an overview of pair
potentials.  Click on the style itself for a full description:</p>
<table border="1" class="docutils">
<colgroup>
<col width="18%" />
<col width="17%" />
<col width="39%" />
<col width="25%" />
</colgroup>
<tbody valign="top">
<tr class="row-odd"><td><a class="reference internal" href="pair_gran.html"><em>bubble</em></a></td>
<td><a class="reference internal" href="pair_gran.html"><em>gran</em></a></td>
<td><a class="reference internal" href="pair_hybrid.html"><em>hybrid</em></a></td>
<td><a class="reference internal" href="pair_hybrid.html"><em>hybrid/overlay</em></a></td>
</tr>
<tr class="row-even"><td><a class="reference internal" href="pair_none.html"><em>none</em></a></td>
<td><a class="reference internal" href="pair_soft.html"><em>soft</em></a></td>
<td><a class="reference internal" href="pair_sph_artvisc_tenscorr.html"><em>sph/artVisc/tensCorr</em></a></td>
<td>&nbsp;</td>
</tr>
</tbody>
</table>
</div>
</div>"""
    lines = lines.split("\n")
    return lines

def buildTheRest():
    lines = """

           </div>
           <div class="articleComments">
            
           </div>
          </div>
          <footer>
  
    <div class="rst-footer-buttons" role="navigation" aria-label="footer navigation">
      
        <a href="Section_gran_models.html" class="btn btn-neutral float-right" title="5. Contact models" accesskey="n" rel="next">Next <span class="fa fa-arrow-circle-right"></span></a>
      
      
        <a href="Section_input_script.html" class="btn btn-neutral" title="3. Input Script" accesskey="p" rel="prev"><span class="fa fa-arrow-circle-left"></span> Previous</a>
      
    </div>
  

  <hr/>

  <div role="contentinfo">
    <p>
        &copy; Copyright 2016, DCS Computing GmbH, JKU Linz and Sandia Corporation.

    </p>
  </div>
  Built with <a href="http://sphinx-doc.org/">Sphinx</a> using a <a href="https://github.com/snide/sphinx_rtd_theme">theme</a> provided by <a href="https://readthedocs.org">Read the Docs</a>. 

</footer>

        </div>
      </div>

    </section>

  </div>
  


  

    <script type="text/javascript">
        var DOCUMENTATION_OPTIONS = {
            URL_ROOT:'./',
            VERSION:'v3.X',
            LANGUAGE:'None',
            COLLAPSE_INDEX:false,
            FILE_SUFFIX:'.html',
            HAS_SOURCE:  true,
            SOURCELINK_SUFFIX: ''
        };
    </script>
      <script type="text/javascript" src="_static/jquery.js"></script>
      <script type="text/javascript" src="_static/underscore.js"></script>
      <script type="text/javascript" src="_static/doctools.js"></script>

  

  
  
    <script type="text/javascript" src="_static/js/theme.js"></script>
  

  
  
  <script type="text/javascript">
      jQuery(function () {
          SphinxRtdTheme.StickyNav.enable();
      });
  </script>
   

</body>
</html>"""
    lines = lines.split("\n")
    return lines

def main():
    lines = buildCommands()
    lines += buildBonds()
    lines += buildComputes()
    lines += buildDumps()
    lines += buildFixes()
    lines += buildPairStyles()
    lines += buildTheRest()

    with open("Section_commands.html",'w') as f:
        for line in lines:
            f.write(line+'\n')
    for line in lines:
        print(line)

if __name__ == "__main__":
    main()