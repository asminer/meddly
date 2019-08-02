
Expectations for developers
----------------------------------------------------------------------

(1) Source code should not contain TAB characters.

(2) There is no mandated coding or naming style, but keep in mind
    that other developers will occasionally need to maintain your
    code, so please keep your code readable.  Use doxygen-style
    comments (similar to Javadoc) for classes and methods.

(3) Ideally, the library builds with NO COMPILER WARNINGs, even with
    -Wall.  This is complicated by the fact that different compilers
    enable different warnings with -Wall, and some code beyond our
    control (other libraries, headers, etc.) may produce warnings.  
    Nonetheless, please work to keep your code free of compiler 
    warnings.

(4) When implementing features or other disruptive changes, the 
    recommended mechanism is to create a branch.  Keep in mind that,
    ideally, the main branch always compiles and works as expected.

(5) Make sure your code passes regression tests (make check) before 
    checking into the main branch

(6) You are strongly encouraged to create your own regression tests
    (that complete in, say, tens of seconds) that run under "make
    check", to test any features you implement.  Otherwise, there is
    a risk that a change to the code will break one of your features.

