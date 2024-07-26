    #version 330
    uniform mat4 MVP;
    layout (location = 0) in vec3 vCol;
    layout (location = 1) in vec3 vPos;
    layout (location = 2) in vec2 aTexCoord;
    layout (location = 3) in vec3 vNormal;
    out vec2 texCoord;
    out vec3 color;
    out vec3 norm;
    out vec3 viewDir;
    uniform vec4 myColor;
    
    uniform mat4 model;
    uniform mat4 view;
    uniform mat4 projection;
    uniform mat4 rotation;

    void main()
    {
        gl_Position = projection * view * model * vec4(vPos, 1.0);
        vec4 wsPos = view * model * vec4(vPos, 1.0);
        color = vCol;
        texCoord = aTexCoord;

        norm = (rotation * vec4(vNormal,1)).xyz;
        vec3 cameraPosition = -view[3].xyz;
        vec3 vertexPositionWorld = (model * vec4(vPos, 1.0)).xyz;
        viewDir = normalize(cameraPosition-vertexPositionWorld);
    };