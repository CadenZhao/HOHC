export default {
  // Cron keep-alive: visit the legacy shinyapps.io app so it is not
  // deleted for inactivity, and check the new site is up.
  async scheduled(event, env, ctx) {
    console.log(JSON.stringify(await ping()));
  },
  // Opening the worker URL runs the same pings manually.
  async fetch(request, env, ctx) {
    return new Response(JSON.stringify(await ping(), null, 2), {
      headers: { "content-type": "application/json" }
    });
  }
};

async function ping() {
  const urls = [
    "https://omicsexplorer.shinyapps.io/HOHC/",
    "https://hohc.pages.dev/"
  ];
  const results = {};
  for (const url of urls) {
    try {
      const r = await fetch(url, { redirect: "follow" });
      await r.text();
      results[url] = r.status;
    } catch (e) {
      results[url] = "ERROR: " + e.message;
    }
  }
  return results;
}
