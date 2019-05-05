#include <stdio.h>
#include <string.h>
#include <sys/time.h>
#include <sys/types.h>
#include <unistd.h>
#include <termios.h>
#include <unistd.h>

static struct termios ori_attr, cur_attr;

int tty_reset(void)
{
    if (tcsetattr(STDIN_FILENO, TCSANOW, &ori_attr) != 0)
    return -1;
    return 0;
}

int tty_set(void)
{

    if ( tcgetattr(STDIN_FILENO, &ori_attr) )
        return -1;

    memcpy(&cur_attr, &ori_attr, sizeof(cur_attr) );
    cur_attr.c_lflag &= ~ICANON;
    // cur_attr.c_lflag |= ECHO;
    cur_attr.c_lflag &= ~ECHO;
    cur_attr.c_cc[VMIN] = 1;
    cur_attr.c_cc[VTIME] = 0;
    if (tcsetattr(STDIN_FILENO, TCSANOW, &cur_attr) != 0)
        return -1;

    return 0;
}

int kbhit(void)
{
    fd_set rfds;
    struct timeval tv;
    int retval;

    /* Watch stdin (fd 0) to see when it has input. */
    FD_ZERO(&rfds);
    FD_SET(0, &rfds);

    /* Wait up to five seconds. */
    tv.tv_sec = 0;
    tv.tv_usec = 0;
    retval = select(1, &rfds, NULL, NULL, &tv);

    /* Don't rely on the value of tv now! */
    if (retval == -1)
    {
        perror("select()");
        return 0;
    }
    else if (retval)
        return 1;

    /* FD_ISSET(0, &rfds) will be true. */
    else
        return 0;

    return 0;
}

/*
int main()
{
    int tty_set_flag;
    tty_set_flag = tty_set();
    while(1)
    {
        if( kbhit() )
        {
            const int key = getchar();
            printf("%c pressed\n", key);
            if(key == 'q')
            break;
        }
        else
        {
        //    fprintf(stderr, "<no key detected>\n");
        }
    }

    if(tty_set_flag == 0)
        tty_reset();

    return 0;
}
*/


