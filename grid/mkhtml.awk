#create an HTML table displaying the CMD images for a particular region
#requires the "pre" argument: two-letter subregion code
BEGIN {
    outfile = pre ".html";
    type[1] = ".dat";
    type[2] = ".mod";
    type[3] = ".chi";
    type[4] = ".dif";
    str[1]  = "Data";
    str[2]  = "Model";
    str[3]  = "Chi";
    str[4]  = "Data-Model";

    ncmd = split( cmds, cmd, "." );
    if ( ncmd==0 ) {
      print "Error: no valid CMD string found.";
      exit(1);
    }
    
    print "<html><head><title>Region " pre "</title></head>" >> outfile;
    print "<body bgcolor=\"\#FFFFFF\">" >> outfile;
    print "<h2>Results for region " pre "</h2>" >> outfile;

    print "<table align=\"center\">" >> outfile;
    print "<tr><td>&nbsp;</td>" >> outfile;

    for ( i=1; i<=ncmd; i++ ) {
      printf("<th>CMD %d</th>", i ) >> outfile;
    }
    printf("\n") >> outfile;

    for ( i=1; i<=4; ++i ) {
	print "<tr>" >> outfile;
	print "<th>" str[i] "</th>" >> outfile;

	for ( j=1; j<=ncmd; ++j ) {
	    print "<td><img src=\"images/" pre type[i] "." cmd[j] ".png\"></td>" >> outfile;
	}
	print "</tr>" >> outfile;
    }

    print "</table>" >> outfile;
    print "</body></html>" >> outfile;
    
}
