    #version 330
    in vec3 color;
    out vec4 fragment;
    uniform vec4 myColor;

    in vec2 texCoord;
    in vec3 norm;
    in vec3 viewDir;
    uniform sampler2D ourTexture;
    uniform sampler2D ourTexture2;
    uniform sampler2D ourTexture3;

    void main()
    {
        float ambientLight = 0.15f;
        float lightBrightNess = 0.8f;
        vec3 lightDir = vec3(0.0,1,0);
        float lightDot = clamp(dot(norm, lightDir), 0.0, 1.0);
        // fragment = vec4(color, 1.0) - myColor;
        // fragment = vec4(1.0,1.0,1.0, 1.0) ;
        vec4 fragmentc = vec4(color, 1.0) - myColor;
        vec4 tex = texture(ourTexture, texCoord*3) * myColor;
        fragment = tex*lightDot*lightBrightNess;
        fragment += tex*ambientLight;
        vec3 halfVector = normalize(lightDir + viewDir);
        float lightDot2 = clamp(dot(norm, halfVector), 0.0, 1.0);
        lightDot2 = pow(lightDot2, texture(ourTexture2, texCoord*3).x * 100);
        //fragment = vec4(lightDot2,lightDot2,lightDot2,1);
        //fragment = vec4(viewDir, 1.0);
        fragment += vec4(lightDot2, lightDot2, lightDot2, 1.0)*lightBrightNess *0.5;
        // fragment = vec4(norm,1);
        // fragment = texture(ourTexture3, texCoord)  ;
    };
