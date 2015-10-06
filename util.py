def wrap(string,width):
    if (width == 0): return string;
    out = "";
    while(string.__len__() > width):
        out += string[:width]+"\n"
        string = string[width:]
    return out
