{
  # run in the source directory where the files are to expanded
  # use recursion to expand include files within include files
  test_include($0)
}

function test_include(str, incl_str, incl_file)
{
  if($1 == "INCLUDE" || $1 == "include")
  {
    incl_str = "1"
    print "C ----------- ^^^^^^^^^^^^ ------------"
    print "C"str
    split(str,infile,"'")
    incl_file = infile[2]
    while((getline < incl_file) > 0)
    {
      test_include($0)
    }
  }
  else
  {
    print str
  }  
  if(incl_str == "1")
  {
    print "C ^^^^^^^^^^^ ------------ ^^^^^^^^^^^"
    close(incl_file)
  }
}
