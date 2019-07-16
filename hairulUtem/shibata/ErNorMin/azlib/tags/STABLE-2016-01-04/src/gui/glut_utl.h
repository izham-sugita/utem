/*********************************************************************
 * glut_utl.h
 *
 * Copyright (C) 2001 AzLib Developers Group
 *
 * Written by
 *  <Kenzen TAKEUCHI> <Mitsunori HIROSE> <Ryu SASAOKA>
 *
 * Refer to log documents for details.
 *********************************************************************/

/* $Id: glut_utl.h 864 2007-01-17 10:02:04Z sasaoka $ */

#ifndef GLUT_UTL_H
#define GLUT_UTL_H

/* キーボード操作用 */
#define KEYBOARD_CTRL_A (1)
#define KEYBOARD_CTRL_B (2)
#define KEYBOARD_CTRL_C (3)
#define KEYBOARD_CTRL_D (4)
#define KEYBOARD_CTRL_E (5)
#define KEYBOARD_CTRL_F (6)
#define KEYBOARD_CTRL_G (7)
#define KEYBOARD_CTRL_H (8)
#define KEYBOARD_CTRL_I (9)
#define KEYBOARD_CTRL_J (10)
#define KEYBOARD_CTRL_K (11)
#define KEYBOARD_CTRL_L (12)
#define KEYBOARD_CTRL_M (13)
#define KEYBOARD_CTRL_N (14)
#define KEYBOARD_CTRL_O (15)
#define KEYBOARD_CTRL_P (16)
#define KEYBOARD_CTRL_Q (17)
#define KEYBOARD_CTRL_R (18)
#define KEYBOARD_CTRL_S (19)
#define KEYBOARD_CTRL_T (20)
#define KEYBOARD_CTRL_U (21)
#define KEYBOARD_CTRL_V (22)
#define KEYBOARD_CTRL_W (23)
#define KEYBOARD_CTRL_X (24)
#define KEYBOARD_CTRL_Y (25)
#define KEYBOARD_CTRL_Z (26)
#define KEYBOARD_BS (8)
#define KEYBOARD_TAB (9)
#define KEYBOARD_RETURN (13)
#define KEYBOARD_ESC (27)
#define KEYBOARD_SPACE (32)
#define KEYBOARD_0 (48)
#define KEYBOARD_1 (49)
#define KEYBOARD_2 (50)
#define KEYBOARD_3 (51)
#define KEYBOARD_4 (52)
#define KEYBOARD_5 (53)
#define KEYBOARD_6 (54)
#define KEYBOARD_7 (55)
#define KEYBOARD_8 (56)
#define KEYBOARD_9 (57)
#define KEYBOARD__A (65)
#define KEYBOARD__B (66)
#define KEYBOARD__C (67)
#define KEYBOARD__D (68)
#define KEYBOARD__E (69)
#define KEYBOARD__F (70)
#define KEYBOARD__G (71)
#define KEYBOARD__H (72)
#define KEYBOARD__I (73)
#define KEYBOARD__J (74)
#define KEYBOARD__K (75)
#define KEYBOARD__L (76)
#define KEYBOARD__M (77)
#define KEYBOARD__N (78)
#define KEYBOARD__O (79)
#define KEYBOARD__P (80)
#define KEYBOARD__Q (81)
#define KEYBOARD__R (82)
#define KEYBOARD__S (83)
#define KEYBOARD__T (84)
#define KEYBOARD__U (85)
#define KEYBOARD__V (86)
#define KEYBOARD__W (87)
#define KEYBOARD__X (88)
#define KEYBOARD__Y (89)
#define KEYBOARD__Z (90)
#define KEYBOARD_A (97)
#define KEYBOARD_B (98)
#define KEYBOARD_C (99)
#define KEYBOARD_D (100)
#define KEYBOARD_E (101)
#define KEYBOARD_F (102)
#define KEYBOARD_G (103)
#define KEYBOARD_H (104)
#define KEYBOARD_I (105)
#define KEYBOARD_J (106)
#define KEYBOARD_K (107)
#define KEYBOARD_L (108)
#define KEYBOARD_M (109)
#define KEYBOARD_N (110)
#define KEYBOARD_O (111)
#define KEYBOARD_P (112)
#define KEYBOARD_Q (113)
#define KEYBOARD_R (104)
#define KEYBOARD_S (115)
#define KEYBOARD_T (116)
#define KEYBOARD_U (117)
#define KEYBOARD_V (118)
#define KEYBOARD_W (119)
#define KEYBOARD_X (120)
#define KEYBOARD_Y (121)
#define KEYBOARD_Z (122)

/* 表示ウインドウ設定 */
#define WINDOW_WIDTH (500)
#define WINDOW_HEIGHT (500)
#define WINDOW_POSITION_X (700)
#define WINDOW_POSITION_Y (50)

/* マウス処理用 */
#define STATE_DOWN (0)
#define STATE_UP   (1)
#define POINT 0
#define POLYGON 1
#define LINE 2
#define MATERIAL_WHITE 0
#define MATERIAL_RED 1
#define MATERIAL_GREEN 2
#define MATERIAL_BLUE 3
#define MATERIAL_CYAN 4
#define MATERIAL_MAGENTA 5
#define MATERIAL_YELLOW 6
#define RESET 98 
#define QUIT 99

#define R(r) (r/255.0) /* Red */
#define G(g) (g/255.0) /* Green */
#define B(b) (b/255.0) /* Blue */
#define H(h) (h/128.0) /* Hi-light */

/* declaration of functions */
RC show_fem_model(ELEMENT_ARRAY element, NODE_ARRAY node);



/* キー配置ナンバー 
key = CTRL+a  Num = 1
key = CTRL+b  Num = 2
key = CTRL+c  Num = 3
key = CTRL+d  Num = 4
key = CTRL+e  Num = 5
key = CTRL+f  Num = 6
key = CTRL+g  Num = 7
key = (BS) CTRL+h     Num = 8 
key = (TAB) CTRL+i    Num = 9
key = CTRL+j  Num = 10
key = CTRL+k  Num = 11
key = CTRL+l  Num = 12
key = (RETURN) CTRL+m Num = 13
key = CTRL+n  Num = 14
key = CTRL+o  Num = 15
key = CTRL+p  Num = 16
key = CTRL+q  Num = 17
key = CTRL+r  Num = 18
key = CTRL+s  Num = 19
key = CTRL+t  Num = 20
key = CTRL+u  Num = 22
key = CTRL+v  Num = 23
key = CTRL+w  Num = 24
key = CTRL+x  Num = 25
key = CTRL+y  Num = 26
key = (ESC) CTRL+z CTRL+[ CTRL+{ Num = 27
key = CTRL+\    Num = 28
key = CTRL+]    Num = 29
key = CTRL+~    Num = 30
key = CTRL+/    Num = 31
key = (SPACE) Num = 32
key = !       Num = 33
key = "       Num = 34
key = #       Num = 35
key = $       Num = 36
key = %       Num = 37
key = &       Num = 38
key = '       Num = 39
key = (       Num = 40
key = )       Num = 41
key = *       Num = 42
key = +       Num = 43
key = ,       Num = 44
key = -       Num = 45
key = .       Num = 46
key = /       Num = 47
key = 0       Num = 48
key = 1       Num = 49
key = 2       Num = 50
key = 3       Num = 51
key = 4       Num = 52
key = 5       Num = 53
key = 6       Num = 54
key = 7       Num = 55
key = 8       Num = 56
key = 9       Num = 57
key = :       Num = 58
key = ;       Num = 59
key = <       Num = 60
key = =       Num = 61
key = >       Num = 62
key = ?       Num = 63
key = @       Num = 64
key = A       Num = 65
key = B       Num = 66
key = C       Num = 67
key = D       Num = 68
key = E       Num = 69
key = F       Num = 70
key = G       Num = 71
key = H       Num = 72
key = I       Num = 73
key = J       Num = 74
key = K       Num = 75
key = L       Num = 76
key = M       Num = 77
key = N       Num = 78
key = O       Num = 79
key = P       Num = 80
key = Q       Num = 81
key = R       Num = 82
key = S       Num = 83
key = T       Num = 84
key = U       Num = 85
key = V       Num = 86
key = W       Num = 87
key = X       Num = 88
key = Y       Num = 89
key = Z       Num = 90
key = [       Num = 91
key = \       Num = 92
key = ]       Num = 93
key = ^       Num = 94
key = _       Num = 95
key = `       Num = 96
key = a       Num = 97
key = b       Num = 98
key = c       Num = 99
key = d       Num = 100
key = e       Num = 101
key = f       Num = 102
key = g       Num = 103
key = h       Num = 104
key = i       Num = 105
key = j       Num = 106
key = k       Num = 107
key = l       Num = 108
key = m       Num = 109
key = n       Num = 110
key = o       Num = 111
key = p       Num = 112
key = q       Num = 113
key = r       Num = 114
key = s       Num = 115
key = t       Num = 116
key = u       Num = 117
key = v       Num = 118
key = w       Num = 119
key = x       Num = 120
key = y       Num = 121
key = z       Num = 122
key = {       Num = 123
key = |       Num = 124
key = }       Num = 125
key = ~       Num = 126
key = (DEL)   Num = 127
*/

#endif /* GLUT_UTL_H */


