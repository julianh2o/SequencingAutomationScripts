def wrap(string,width):
    out = "";
    while(string.__len__() > width):
        out += string[:width]+"\n"
        string = string[width:]
    return out
