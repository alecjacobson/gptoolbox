These files patch MATLAB to automatically remember ALL commands in your history
(rather than just the last 20KB's worth)

On Unix, issue:
(You may have to change R2014a to whatever matlab version you have)

    cp ~/.matlab/R2014a/history.m ~/.matlab/R2014a/history_saved.m
    mv preserve_history.m ~/.matlab/R2014a/preserve_history.m

If you don't already have a `finish.m` or `startup.m` you can issue:

    mv finish.m ~/.matlab/R2014a/finish.m
    mv startup.m ~/.matlab/R2014a/startup.m

Otherwise you'll have to append these files to your existing `finish.m` and
`startup.m` files.



