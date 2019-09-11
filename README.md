# The Planet

```
// GROUP NUMBER:13
// STUDENT NAME:WU PENGYU
// NUS User ID.: t0918566
// STUDENT NAME: WANG SIYUAN
// NUS User ID.: t0918594
// STUDENT NAME: FAN QIANYI	
// NUS User ID.: t0918683
```
## VIEW


We create a shader to render a solar system. In this shader, we implement a planet class with rusty surface which phone shading supported. We also implement a method to render the trail of each planet.

If you have configure the enviroment right, you will see following scene.

At first, you can see the whole solar system. You can use the direction keys ⬆️⤵️⬅️➡️ to move the position of camera.

![db139cbef1e9921a6492193121e17777](The Planet.resources/DB139CBEF1E9921A6492193121E17777.jpg)


After twenty or thirty seconds, the scene will shrink to the sun. You enter a world named Mirror Space. Well, it's just an animation。

![05a29df62d8c5898fc9fe8ef0b431004](The Planet.resources/05A29DF62D8C5898FC9FE8EF0B431004.jpg)


That's all.

## CONFIG

Here is the configure procedure.

We create three buffer in shadertoy, named buffer A, buffer B and Image. Each buffer have its own texture configuration.

After pasting the each codes from files to the corresponding buffer, we need to configurate the texture it reads from.
In buffer A, you should select which texture to read in each iChannel. You may configurate iChannel like this.

![b76a557cc96c720c8273578c4d2f4bf3](The Planet.resources/B76A557CC96C720C8273578C4D2F4BF3.png)

In buffer B, you needn't change any thing. Just left all the iChannel blank.

In the buffer of image, You may configurate iChannel like this.

![f0d20e13144301b5f3081df2e11be3c0](The Planet.resources/F0D20E13144301B5F3081DF2E11BE3C0.jpg)

So, finally the environment looks like this:

![1bcaf7cb3819f01c30234dbd226de430](The Planet.resources/屏幕快照 2019-07-24 下午11.47.22.png)

![52aef29142d1a812013face1a17fb96b](The Planet.resources/屏幕快照 2019-07-24 下午11.47.26.png)

![577d6562657cacbda4012dc825011c73](The Planet.resources/屏幕快照 2019-07-24 下午11.47.30.png)

If it doesn't work well, Here is the temp link for the shader.

https://www.shadertoy.com/view/tllSzs