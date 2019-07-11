<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
    <xsl:output omit-xml-declaration="yes" version="1.0" indent="no" method="text" encoding="UTF-8" media-type="text/plain" />
    <xsl:strip-space elements="*" />
    <xsl:preserve-space elements="para"/>

    <!-- read all doxygenindex/compound entries in the index file and process referenced files individually -->
    <xsl:template match="/">
        <xsl:text># SIM5 Library Reference</xsl:text>
        <xsl:text>&#xa;&#xa;</xsl:text>
        <xsl:text>SIM5 is a collection of C routines for relativistic raytracing and radiation transfer. It has a special focus on raytracing from accretion disks, tori, hot spots or any custom 3D configuration of matter in Kerr geometry, but it can be used with any other metrics as well. It can handle both optically thick and thin sources as well as transport of polarization properties and helps to calculate the propagation of light rays from the source to an observer through a curved spacetimes. The library is threas-safe (with a few documented exceptions) compiles with Nvidia CUDA compiler which opens the door to massive parallelization using GPUs.</xsl:text>
        <xsl:text>&#xa;&#xa;</xsl:text>
        <xsl:text>SIM5 also comes with a Python interface making it very easy to call its functions directly from Python scripts. In addition it provides few Python classes to handle some more complex tasks.</xsl:text>
                 
        <xsl:text>&#xa;&#xa;</xsl:text>
        <xsl:text>The following documentation provides a detailed reference to the functions of the library. The library also comes with couple of examples that illustrate how to piece the individual routines together to do some useful stuff. Routines are grouped together according to specific topics. The content gives a list of these topics with link to corresponding parts of the docs. If you are looking for something specific, use Ctrl-F.</xsl:text>
                 
        <xsl:text>&#xa;&#xa;</xsl:text>
        <xsl:text><![CDATA[<br>]]>&#xa;&#xa;</xsl:text>

        <xsl:text>## Content&#xa;</xsl:text>
        <xsl:for-each select="/doxygenindex/compound[@kind='file']">
            <xsl:variable name="filename" select="substring-before(document(concat('xml/',@refid,'.xml'))/doxygen/compounddef/compoundname, '.')" />
            <xsl:variable name="filebrief" select="substring-before(document(concat('xml/',@refid,'.xml'))/doxygen/compounddef/briefdescription/para, '.')" />
            <xsl:text>* [</xsl:text><xsl:value-of select="$filebrief"/><xsl:text> (</xsl:text><xsl:value-of select="$filename"/><xsl:text>)](#</xsl:text><xsl:value-of select="$filename"/><xsl:text>)&#xa;</xsl:text>
        </xsl:for-each>
        <xsl:text>&#xa;&#xa;</xsl:text>
        <xsl:text><![CDATA[<br>]]>&#xa;&#xa;</xsl:text>

        <xsl:for-each select="/doxygenindex/compound[@kind='file']">
            <xsl:apply-templates select="document(concat('xml/',@refid,'.xml'))/doxygen/compounddef" />
        </xsl:for-each>

    <!--

Content
* [sim5kerr.c][]
* [sim5disk-nt.c][]
        </xsl:text>
        -->
        <!--<xsl:apply-templates select="/doxygenindex/compound" />-->
        <!--<xsl:apply-templates select="document('xml/sim5disk-nt_8c.xml')/doxygen/compounddef"/>-->
    </xsl:template>
    



    <!-- read all doxygenindex/compound entries in the index file and process referenced files individually -->
    <!--
    <xsl:template match="/doxygenindex/compound">
        <xsl:apply-templates select="document('xml/sim5disk-nt_8c.xml')/doxygen/compounddef"/>
    </xsl:template>
    -->



    <!-- process definitions for of each file section -->
    <xsl:template match="compounddef[@kind='file']">
        <!-- file description -->
        <xsl:text><![CDATA[<br>]]>&#xa;&#xa;</xsl:text>
        <xsl:text>## </xsl:text>
        <xsl:text><![CDATA[<a name="]]></xsl:text><xsl:value-of select="substring-before(compoundname, '.')" /><xsl:text><![CDATA["></a>]]></xsl:text>
        <xsl:value-of select="compoundname" />
        <xsl:text> - </xsl:text>

        <xsl:apply-templates select="briefdescription">
            <xsl:with-param name='make-para' select="0"/>
            <xsl:with-param name='remove-dot' select="1"/>
        </xsl:apply-templates>
        <xsl:text>&#xa;&#xa;</xsl:text>

        <!-- detailed description -->
        <xsl:apply-templates select="detaileddescription"/>
        <xsl:text>&#xa;</xsl:text>

        <!-- process all file sections -->
        <xsl:apply-templates select="sectiondef" />

        <!--
        <xsl:apply-templates select="listofallmembers"/>
        Location: [[\\<xsl:value-of select="./location/@file"/>]]
        -->
        <xsl:text><![CDATA[<br><br><br>]]>&#xa;&#xa;</xsl:text>
    </xsl:template>

    <xsl:template match="compounddef">
    </xsl:template>

<!--

    <xsl:template match="sectiondef[@kind='func']">
        <!- -### Functions<xsl:text>&#xa;</xsl:text> - ->
    </xsl:template>
-->    

    <xsl:template match="sectiondef[@kind='func']">
        <!--<xsl:text>### Functions&#xa;</xsl:text>-->
        <xsl:apply-templates select="memberdef" />
    </xsl:template>

    <xsl:template match="sectiondef[@kind='var']">
        <!--
        <xsl:text>### Variables&#xa;</xsl:text>
        <xsl:text>vars follow...&#xa;</xsl:text>
        -->
    </xsl:template>

    <xsl:template match="sectiondef[@kind='define']">
        <!--<xsl:text>### Defines&#xa;</xsl:text>-->

        <xsl:text>| Name | Description | Value |&#xa;</xsl:text>
        <xsl:text>|------|-------------|-------|&#xa;</xsl:text>
        <xsl:for-each select="memberdef[@kind='define']">
            <xsl:text>| </xsl:text>
            <xsl:value-of select="name" />
            <xsl:text> | </xsl:text>
            <xsl:apply-templates select="briefdescription">
                <xsl:with-param name='print-as-para' select="0"/>
                <xsl:with-param name='remove-dot' select="1"/>
            </xsl:apply-templates>

            <xsl:if test="detaileddescription/para">
            <xsl:text><![CDATA[<br>]]></xsl:text>
                <xsl:apply-templates select="detaileddescription/para" />
            </xsl:if>
        
            <xsl:text> | </xsl:text>
            <xsl:value-of select="initializer" />
            <xsl:text> |&#xa;</xsl:text>
        </xsl:for-each>


    </xsl:template>



    <!--
    <xsl:template match="sectiondef[@kind  = 'private-attrib']">
        ## Private Attributes
        <xsl:apply-templates select="memberdef" />
    </xsl:template>

    <xsl:template match="sectiondef[@kind  = 'property']">
        ## Properties
        <xsl:apply-templates select="memberdef" />
    </xsl:template>

    <xsl:template match="sectiondef[@kind  = 'private-static-func']">
        ## Private Static Functions
        <xsl:apply-templates select="memberdef" />
    </xsl:template>


    <xsl:template match="listofallmembers">
        ## All Members
        <xsl:apply-templates select="member" />
    </xsl:template>


    <xsl:template match="inheritancegraph">
        |
        <xsl:for-each select="./node">
            [[#
            <xsl:value-of select="./label" />
            ]]
            <xsl:choose>
                <xsl:when test="position() != last()">-&gt;</xsl:when>
            </xsl:choose>
        </xsl:for-each>
        |
    </xsl:template>
    -->

    <xsl:template match="memberdef[@kind='function']">
        <xsl:text>#### </xsl:text><xsl:value-of select="./name" /><xsl:text>()&#xa;&#xa;</xsl:text>
        <xsl:apply-templates select="briefdescription">
            <xsl:with-param name='make-para' select="1"/>
            <xsl:with-param name='remove-dot' select="0"/>
        </xsl:apply-templates>
        <xsl:text>&#xa;</xsl:text>
        <xsl:text>    </xsl:text><xsl:value-of select="type" /><xsl:text> </xsl:text><xsl:value-of select="name" /><xsl:value-of select="argsstring" /><xsl:text>&#xa;</xsl:text>
        <xsl:text>&#xa;</xsl:text>
        <xsl:apply-templates select="detaileddescription" />
        <xsl:text><![CDATA[<br>]]>&#xa;&#xa;</xsl:text>
        <xsl:text>---------&#xa;&#xa;</xsl:text>
    </xsl:template>


    <xsl:template match="memberdef">
    </xsl:template>

    

<!--
    <xsl:template match="memberdef[(@kind = 'variable') or (@kind = 'property')]">
        ^
        <xsl:value-of select="./type" />
        ^**[[#
        <xsl:value-of select="./name" />
        ]]**^
|//
        <xsl:value-of select="./briefdescription" />
        // ||
    </xsl:template>


    <xsl:template match="member">
        <xsl:call-template name="detaileddescription">
            <xsl:with-param name="nodes" select="../../sectiondef/memberdef[@id=current()/@refid]" />
        </xsl:call-template>
    </xsl:template>


    <xsl:template name="detaileddescription">
        <xsl:param name="nodes" />
        <xsl:for-each select="$nodes">
            ====
            <xsl:value-of select="./name" />
            ==== 
^
            <xsl:value-of select="./definition" />
            ^
            <xsl:value-of select="./argsstring" />
            ^
| //
            <xsl:value-of select="./briefdescription" />
            // ||
            <xsl:apply-templates select="./detaileddescription" />
        </xsl:for-each>
    </xsl:template>
-->


    <!-- process memberdef/briefdescription paragraph in boldface 
    <xsl:template match="memberdef/briefdescription">
        <!- -<xsl:text>**</xsl:text>- ->
        <xsl:apply-templates select="."/>
        <!- -<xsl:text>**</xsl:text>- ->
        <xsl:text>&#xa;&#xa;</xsl:text>
    </xsl:template>
    -->


    <!-- process briefdescription -->
    <xsl:template match="briefdescription">
        <!--
        there are two params to the template:
          - remove-dot: controls whether/not to remove trailing dot character ('.') [default=YES]
          - print-as-para: controls whether/not to print the text as a separate paragraph [default=YES]
        -->
        <xsl:param name="remove-dot" select="1"/>
        <xsl:param name="print-as-para" select="1"/>

        <!-- get text string and its length -->
        <xsl:variable name="text" select="normalize-space(para)"/>
        <xsl:variable name="length" select="string-length($text)"/>

        <xsl:if test="$length > 0">
            <xsl:choose>
                <xsl:when test="$remove-dot">
                    <!-- print without the trailing dot -->
                    <xsl:value-of select="substring($text,1,$length - 1)" />
                    <xsl:if test="substring($text,$length) != '.'">
                        <xsl:value-of select="substring($text,$length,1)" />
                    </xsl:if>
                </xsl:when>
                <xsl:otherwise>
                    <!-- print with the trailing dot -->
                    <xsl:value-of select="$text" />
                    <xsl:if test="substring($text,$length) != '.'">
                        <xsl:text>.</xsl:text>
                    </xsl:if>
                </xsl:otherwise>
            </xsl:choose>
            
            <xsl:if test="$print-as-para">
                <xsl:text>&#xa;&#xa;</xsl:text>
            </xsl:if>
        </xsl:if>
    </xsl:template>
    
    
    <xsl:template match="detaileddescription">
        <xsl:for-each select="para">
            <xsl:apply-templates/>
            <xsl:text>&#xa;</xsl:text>
            <xsl:text>&#xa;</xsl:text>
        </xsl:for-each>
    </xsl:template>


    <xsl:template match="para">
        <xsl:param name="level"/>
        <xsl:apply-templates>
            <xsl:with-param name="level" select="$level"/>
        </xsl:apply-templates>
    </xsl:template>


    <xsl:template match="detaileddescription/para/parameterlist[@kind='param']">
        <xsl:text>**Parameters**</xsl:text><xsl:text>&#xa;&#xa;</xsl:text>
        <xsl:for-each select="parameteritem">
            <xsl:text>* </xsl:text>
            <xsl:text>**</xsl:text><xsl:value-of select="parameternamelist/parametername"/><xsl:text>**</xsl:text>
            <xsl:text>: </xsl:text>
            <xsl:apply-templates select="parameterdescription/para">
                <xsl:with-param name="level" select="2"/>
            </xsl:apply-templates>
            <xsl:text>&#xa;</xsl:text>
        </xsl:for-each>
        <xsl:text>&#xa;</xsl:text>
    </xsl:template>


    <xsl:template match="detaileddescription/para/simplesect[@kind='return']">
        <xsl:text>**Return value**</xsl:text><xsl:text>&#xa;&#xa;</xsl:text>
        <xsl:value-of select="para" /><xsl:text>&#xa;</xsl:text>
        <xsl:text>&#xa;</xsl:text>
    </xsl:template>


    <xsl:template match="formula">
        <xsl:variable name="eq0" select="." />
        <xsl:choose>
            <xsl:when test="contains($eq0,'\[')">
                <!--
                <xsl:variable name="eq1" select="concat(substring-before($eq0,'\['), '$$', substring-after($eq0,'\['))" />
                <xsl:variable name="eq2" select="concat(substring-before($eq1,'\]'), '$$', substring-after($eq1,'\['))" />
                <xsl:text>&#xa;&#xa;</xsl:text><xsl:value-of select="$eq2"/><xsl:text>&#xa;&#xa;</xsl:text>
                -->
                <!-- the following implements GitLab KaTex syantax for stand-alone equations -->
                <xsl:variable name="eq1" select="concat(substring-before($eq0,'\['), substring-after($eq0,'\['))" />
                <xsl:variable name="eq2" select="concat(substring-before($eq1,'\]'), substring-after($eq1,'\['))" />
                <xsl:text>&#xa;```math&#xa;</xsl:text><xsl:value-of select="$eq2"/><xsl:text>&#xa;```&#xa;</xsl:text>
            </xsl:when>
            <xsl:otherwise>
                <!-- the following implements GitLab KaTex syantax for in-line equations -->
                <xsl:variable name="eq1" select="concat(substring-before($eq0,'$'), substring-after($eq0,'$'))" />
                <xsl:variable name="eq2" select="concat(substring-before($eq1,'$'), substring-after($eq1,'$'))" />
                <xsl:text>$`</xsl:text><xsl:value-of select="$eq2"/><xsl:text>`$</xsl:text>
            </xsl:otherwise>
        </xsl:choose>
    </xsl:template>


    <xsl:template match="verbatim">
        <xsl:text>&#xa;</xsl:text>
        <xsl:text>```&#xa;</xsl:text>
        <xsl:value-of select="."/>
        <xsl:text>```&#xa;</xsl:text>
        <xsl:text>&#xa;</xsl:text>
    </xsl:template>


    <xsl:template match="computeroutput">
        <xsl:text>`</xsl:text><xsl:value-of select="."/><xsl:text>`</xsl:text>
    </xsl:template>
    
    
    <xsl:template match="itemizedlist">
        <xsl:param name="level" select="1"/>
        <xsl:text>&#xa;</xsl:text>
        <xsl:for-each select="./listitem">
            <xsl:if test="$level > 1">
                <xsl:text>  </xsl:text>
            </xsl:if>
             <xsl:text>* </xsl:text><xsl:apply-templates/><xsl:text>&#xa;</xsl:text>
        </xsl:for-each>
        <xsl:text>&#xa;</xsl:text>
    </xsl:template>


    <xsl:template name="string-replace-all">
        <xsl:param name="text"/>
        <xsl:param name="replace"/>
        <xsl:param name="by"/>
        <xsl:choose>
            <xsl:when test="contains($text,$replace)">
                <xsl:value-of select="substring-before($text,$replace)"/>
                <xsl:value-of select="$by"/>
                <xsl:call-template name="string-replace-all">
                    <xsl:with-param name="text" select="substring-after($text,$replace)"/>
                    <xsl:with-param name="replace" select="$replace"/>
                    <xsl:with-param name="by" select="$by"/>
                </xsl:call-template>
            </xsl:when>
        <xsl:otherwise>
            <xsl:value-of select="$text"/>
        </xsl:otherwise>
    </xsl:choose>



</xsl:template>



</xsl:stylesheet>
