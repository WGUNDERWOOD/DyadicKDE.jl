using PyPlot

# start plot
plt.ioff()
(fig, ax) = plt.subplots(figsize=(5, 5))
for side in ["bottom", "top", "left", "right"]
    ax.spines[side].set_color("#FFFFFF00")
end
plt.xticks([])
plt.yticks([])

# colors
col_dark = "#111111"
col_light = "#f9f9ee"
col_purple = "#d5b8f5"
col_1 = "#a8f896"
col_2 = "#2e95ef"
col_3 = "#ffec02"
col_4 = "#f56666"

# outer box
lw = 12
x_eps = 0.037
y_eps = 0.020
top = 0.59
xs_outer = [0.5, 1 + x_eps, 1 + x_eps, -x_eps, -x_eps, 0.5]
ys_outer = [-y_eps, -y_eps, top + y_eps, top + y_eps, -y_eps, -y_eps]
outer_col = col_purple
plt.plot(xs_outer, ys_outer, c=outer_col, lw=2 * lw, zorder=0)
plt.fill(xs_outer, ys_outer, c=outer_col, lw=0)

# upper curve
n = 100
xs = range(0, 1, length=n)
ys = [0.6 * (1.2 * x - 0.63)^3 + 0.67 - 0.3 * x
      for x in xs]
all_xs = [[0.5, 0]; xs; [1, 0.5]]
upper_ys = [[0, 0]; ys; [0, 0]]
plt.plot(all_xs, upper_ys, c=col_dark, lw=lw)

# lower curve
n = 100
a = 0.25
ys = [0.3 * (1.2 * (x - a) - 0.63)^3 + 0.48 - 0.3 * (x - a)
      for x in xs]
all_xs = [[0.5, 0]; xs; [1, 0.5]]
lower_ys = [[0, 0]; ys; [0, 0]]
plt.plot(all_xs, lower_ys, c=col_dark, lw=lw)
plt.fill_between(all_xs, lower_ys, upper_ys, fc=col_1)
plt.fill_between(all_xs, 0, lower_ys, fc=col_light)

# network nodes
v1 = [0.18, 0.1]
v2 = [0.35, 0.3]
v3 = [0.81, 0.13]

# plot network edges
plt.plot([v1[1], v2[1]], [v1[2], v2[2]], c=col_dark, lw=7, zorder=10)
plt.plot([v2[1], v3[1]], [v2[2], v3[2]], c=col_dark, lw=7, zorder=10)
plt.plot([v3[1], v1[1]], [v3[2], v1[2]], c=col_dark, lw=7, zorder=10)

# plot network node buffers
r = 4000
ax.scatter([v1[1]], [v1[2]], s=r, lw=0, fc=col_light, zorder=11)
ax.scatter([v2[1]], [v2[2]], s=r, lw=0, fc=col_light, zorder=11)
ax.scatter([v3[1]], [v3[2]], s=r, lw=0, fc=col_light, zorder=11)

# plot network nodes
r = 2000
ax.scatter([v1[1]], [v1[2]], s=r, lw=6, fc=col_2, ec=col_dark, zorder=12)
ax.scatter([v2[1]], [v2[2]], s=r, lw=6, fc=col_3, ec=col_dark, zorder=12)
ax.scatter([v3[1]], [v3[2]], s=r, lw=6, fc=col_4, ec=col_dark, zorder=12)

# save
plt.savefig("../docs/src/assets/logo.svg", dpi=1000, transparent=true)
plt.close("all")
